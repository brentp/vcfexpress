use mlua::Lua;
use rust_htslib::bcf::{
    self,
    header::{TagLength, TagType},
    Read,
};
use std::{collections::HashMap, hash::Hash, io::Write};

use crate::variant::Variant;

pub struct VCFExpr<'lua> {
    lua: &'lua Lua,
    vcf_reader: Option<bcf::Reader>,
    template: Option<mlua::Function<'lua>>,
    writer: Option<EitherWriter>,
    expressions: Vec<mlua::Function<'lua>>,
    set_expressions: HashMap<InfoFormat, ((TagType, TagLength), mlua::Function<'lua>)>,
    globals: mlua::Table<'lua>,
    variants_evaluated: usize,
    variants_passing: usize,
}

pub enum StringOrVariant {
    String(String),
    Variant(Option<bcf::Record>),
    None,
}

pub enum EitherWriter {
    Vcf(bcf::Writer),
    File(std::io::BufWriter<std::fs::File>),
    Stdout(std::io::BufWriter<std::io::Stdout>),
}

impl EitherWriter {
    pub fn translate(&mut self, record: &mut bcf::Record) {
        if let EitherWriter::Vcf(ref mut w) = self {
            w.translate(record);
        }
    }

    pub fn write(&mut self, sob: &mut StringOrVariant) -> std::io::Result<()> {
        match sob {
            StringOrVariant::None => Ok(()),
            StringOrVariant::Variant(None) => Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                "expected VCF record got None",
            )),
            StringOrVariant::Variant(Some(ref mut record)) => {
                if let EitherWriter::Vcf(ref mut wtr) = self {
                    match wtr.write(record) {
                        Ok(_) => Ok(()),
                        Err(e) => Err(std::io::Error::new(std::io::ErrorKind::Other, e)),
                    }
                } else {
                    // error because we should not be writing a record to a file or stdout
                    Err(std::io::Error::new(
                        std::io::ErrorKind::Other,
                        "expected VCF writer without template",
                    ))
                }
            }
            StringOrVariant::String(s) => match self {
                EitherWriter::Vcf(ref mut _wtr) => Err(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    "did not VCF writer with template",
                )),
                EitherWriter::File(ref mut f) => writeln!(f, "{}", s),
                EitherWriter::Stdout(ref mut f) => writeln!(f, "{}", s),
            },
        }
    }
}

fn get_vcf_format(path: &str) -> bcf::Format {
    if path.ends_with(".bcf") || path.ends_with(".bcf.gz") {
        bcf::Format::Bcf
    } else {
        bcf::Format::Vcf
    }
}

fn process_template(template: Option<String>, lua: &Lua) -> Option<mlua::Function<'_>> {
    if let Some(template) = template.as_ref() {
        // check if template contains backticks
        let return_pre = if template.contains("return ") {
            ""
        } else {
            "return "
        };
        // add the backticks and return if needed.
        let expr = if template.contains('`') {
            format!("{}{}", return_pre, template)
        } else {
            format!("{} `{}`", return_pre, template)
        };
        Some(lua.load(expr).into_function().expect("error in template"))
    } else {
        None
    }
}

#[derive(Debug, PartialEq, Eq, Hash)]
enum InfoFormat {
    Info(String),
    #[allow(dead_code)]
    Format(String),
}

#[derive(Debug)]
enum InfoFormatValue {
    Bool(bool),
    Float(f32),
    Integer(i32),
    String(String),
}

impl<'lua> VCFExpr<'lua> {
    /// Create a new VCFExpr object. This object will read a VCF file, evaluate a set of expressions.
    /// The expressions should return a boolean. Evaluations will stop on the first true expression.
    /// If a template is provided, the template will be evaluated in the same scope as the expression and used
    /// to generate the text output. If no template is provided, the VCF record will be written to the output.
    pub fn new(
        lua: &'lua Lua,
        vcf_path: String,
        expression: Vec<String>,
        set_expression: Vec<String>,
        template: Option<String>,
        lua_prelude: Option<String>,
        output: Option<String>,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        lua.load(crate::pprint::PPRINT).set_name("pprint").exec()?;
        lua.load(crate::pprint::PRELUDE)
            .set_name("prelude")
            .exec()?;

        let mut reader = match vcf_path.as_str() {
            "-" | "stdin" => bcf::Reader::from_stdin()?,
            _ => bcf::Reader::from_path(&vcf_path)?,
        };
        _ = reader.set_threads(2);
        crate::register(lua)?;
        let globals = lua.globals();
        let template = process_template(template, lua);

        let exps: Vec<_> = expression
            .iter()
            .map(|exp| {
                lua.load(exp)
                    .set_name(exp)
                    .into_function()
                    .expect("error in expression")
            })
            .collect();

        let mut hv = bcf::header::HeaderView::new(unsafe {
            rust_htslib::htslib::bcf_hdr_dup(reader.header().inner)
        });

        if let Some(lua_code) = lua_prelude {
            let code = std::fs::read_to_string(lua_code)?;
            lua.scope(|scope| {
                globals.raw_set("header", scope.create_any_userdata_ref_mut(&mut hv)?)?;
                lua.load(&code).exec()
            })?;
        }

        let info_exps = VCFExpr::load_info_expressions(lua, &mut hv, set_expression)?;

        let header = bcf::header::Header::from_template(&hv);

        let writer = if template.is_none() {
            EitherWriter::Vcf(if let Some(output) = output {
                let format = get_vcf_format(&output);
                let mut wtr =
                    bcf::Writer::from_path(&output, &header, !output.ends_with(".gz"), format)?;
                _ = wtr.set_threads(2);
                wtr
            } else {
                bcf::Writer::from_stdout(&header, true, bcf::Format::Vcf)?
            })
        } else if output.is_none() || output.as_ref().unwrap() == "-" {
            EitherWriter::Stdout(std::io::BufWriter::new(std::io::stdout()))
        } else {
            let file = std::fs::File::create(output.unwrap())?;
            EitherWriter::File(std::io::BufWriter::new(file))
        };

        Ok(VCFExpr {
            lua,
            vcf_reader: Some(reader),
            template,
            writer: Some(writer),
            expressions: exps,
            set_expressions: info_exps,
            globals,
            variants_evaluated: 0,
            variants_passing: 0,
        })
    }

    #[allow(clippy::type_complexity)]
    fn load_info_expressions(
        lua: &'lua Lua,
        hv: &mut bcf::header::HeaderView,
        info_expressions: Vec<String>,
    ) -> Result<
        HashMap<InfoFormat, ((TagType, TagLength), mlua::Function<'lua>)>,
        Box<dyn std::error::Error>,
    > {
        let info_exps: HashMap<_, _> = info_expressions
            .iter()
            .map(|exp| {
                let name_exp = exp
                    .split_once('=')
                    .expect("invalid info expression should have name=$expression");
                let t = hv
                    .info_type(name_exp.0.as_bytes())
                    .unwrap_or_else(|_| panic!("ERROR: info field '{}' not found. Make sure it was added to the header in prelude if needed.", name_exp.0));
                (
                    InfoFormat::Info(name_exp.0.to_string()),
                    (
                        t,
                        lua.load(name_exp.1)
                            .set_name(exp)
                            .into_function()
                            .unwrap_or_else(|_| panic!("error in expression: {}", exp)),
                    ),
                )
            })
            .collect();
        Ok(info_exps)
    }

    /// Add lua code to the Lua interpreter. This code will be available to the expressions and the template.
    /// These are not the variant expressions, but rather additional Lua code that can be used as a library.
    pub fn add_lua_code(&mut self, path: &str) -> Result<(), Box<dyn std::error::Error>> {
        let code = std::fs::read_to_string(path)?;
        match self.lua.load(&code).set_name(path).exec() {
            Ok(_) => (),
            Err(e) => {
                log::error!("Error loading Lua code from {}: {}", path, e);
                return Err(e.into());
            }
        }
        Ok(())
    }

    /// Return a reference to the bcf::Reader object.
    pub fn reader(&mut self) -> bcf::Reader {
        self.vcf_reader.take().expect("reader already taken")
    }

    pub fn writer(&mut self) -> EitherWriter {
        self.writer.take().expect("writer already taken")
    }

    // this is called from in the scope and lets us evaluate the info expressions.
    // we collect the results to be used outside the scope where we can get a mutable variant.
    fn evaluate_info_expressions(
        &self,
        info_results: &mut HashMap<String, InfoFormatValue>,
    ) -> mlua::Result<()> {
        for (inf, ((tagtyp, _taglen), expr)) in self.set_expressions.iter() {
            if let InfoFormat::Info(tag) = inf {
                let t = match tagtyp {
                    TagType::Flag => {
                        let b = expr.call::<_, bool>(())?;
                        InfoFormatValue::Bool(b)
                    }
                    TagType::Float => {
                        let f = expr.call::<_, f32>(())?;
                        InfoFormatValue::Float(f)
                    }
                    TagType::Integer => {
                        let i = expr.call::<_, i32>(())?;
                        InfoFormatValue::Integer(i)
                    }
                    TagType::String => {
                        let s = expr.call::<_, String>(())?;
                        InfoFormatValue::String(s)
                    }
                };
                info_results.insert(tag.clone(), t);
            }
        }
        Ok(())
    }

    pub fn evaluate(&mut self, record: bcf::Record) -> std::io::Result<StringOrVariant> {
        let mut variant = Variant::new(record);
        self.variants_evaluated += 1;
        let mut info_results = HashMap::new();
        let eval_result = self.lua.scope(|scope| {
            let ud = match scope.create_any_userdata_ref_mut(&mut variant) {
                Ok(ud) => ud,
                Err(e) => return Err(e),
            };
            match self.globals.raw_set("variant", ud) {
                Ok(_) => (),
                Err(e) => return Err(e),
            }
            self.evaluate_info_expressions(&mut info_results)?;
            // we have many expressions, we stop on the first passing expression. The result of this scope
            // can be either a bool, or a string (if we have a template).
            for exp in &self.expressions {
                match exp.call::<_, bool>(()) {
                    Err(e) => return Err(e),
                    Ok(true) => {
                        self.variants_passing += 1;
                        if let Some(template) = &self.template {
                            // if we have a template, we want to evaluate it in this same scope.
                            return match template.call::<_, String>(()) {
                                Ok(res) => Ok(StringOrVariant::String(res)),
                                Err(e) => {
                                    log::error!("Error in template: {}", e);
                                    return Err(e);
                                }
                            };
                        }
                        return Ok(StringOrVariant::Variant(None));
                    }
                    Ok(false) => {}
                }
            }

            Ok(StringOrVariant::None)
        });

        let mut record = variant.take();
        for (stag, value) in info_results {
            let tag = stag.as_bytes();
            //debug!("Setting info field: {}: {:?}", stag, value);
            let result = match value {
                InfoFormatValue::Bool(b) => {
                    if b {
                        record.push_info_flag(tag)
                    } else {
                        record.clear_info_flag(tag)
                    }
                }
                InfoFormatValue::Float(f) => record.push_info_float(b"af_copy", &[f]),
                InfoFormatValue::Integer(i) => record.push_info_integer(tag, &[i]),
                InfoFormatValue::String(s) => record.push_info_string(tag, &[s.as_bytes()]),
            };
            match result {
                Ok(_) => (),
                Err(e) => {
                    log::error!("Error setting info field: {}: {}", stag, e);
                    return Err(std::io::Error::new(std::io::ErrorKind::Other, e));
                }
            }
        }
        match eval_result {
            Ok(StringOrVariant::Variant(None)) => Ok(StringOrVariant::Variant(Some(record))),
            Ok(b) => Ok(b),
            Err(e) => Err(std::io::Error::new(std::io::ErrorKind::Other, e)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use mlua::Lua;

    #[test]
    fn test_process_template_with_none() {
        let lua = Lua::new();
        assert_eq!(process_template(None, &lua), None);
    }

    #[test]
    fn test_process_template_with_backticks() {
        let lua = Lua::new();
        let template = Some("`print('Hello, World!')`".to_string());
        let result = process_template(template, &lua);
        assert!(result.is_some());
    }

    #[test]
    fn test_process_template_without_backticks() {
        let lua = Lua::new();
        let template = Some("print('Hello, World!')".to_string());
        let result = process_template(template, &lua);
        assert!(result.is_some());
        // execute the result
        let result = result.unwrap();
        let result = result.call::<_, String>(());
        assert!(result.is_ok());
    }

    #[test]
    fn test_process_template_with_return() {
        let lua = Lua::new();
        let template = Some("return `42`".to_string());
        let result = process_template(template, &lua);
        assert!(result.is_some());
        let result = result.unwrap();
        let result = result.call::<_, i32>(());
        if let Ok(result) = result {
            assert_eq!(result, 42);
        } else {
            panic!("error in template");
        }
    }

    #[test]
    #[should_panic(expected = "error in template")]
    fn test_process_template_with_invalid_lua() {
        let lua = Lua::new();
        let template = Some("return []invalid_lua_code".to_string());
        process_template(template, &lua);
    }
}
