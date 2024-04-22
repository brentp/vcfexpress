use mlua::Lua;
use rust_htslib::bcf::{self, Read};
use std::io::Write;

use crate::variant::Variant;

pub struct VCFExpr<'lua> {
    lua: &'lua Lua,
    vcf_reader: bcf::Reader,
    template: Option<mlua::Function<'lua>>,
    writer: EitherWriter,
    expressions: Vec<mlua::Function<'lua>>,
    globals: mlua::Table<'lua>,
    variants_evaluated: usize,
    variants_passing: usize,
}

pub enum StringOrBool {
    String(String),
    Bool(bool),
}

enum EitherWriter {
    Vcf(bcf::Writer),
    File(std::io::BufWriter<std::fs::File>),
    Stdout(std::io::BufWriter<std::io::Stdout>),
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

impl<'lua> VCFExpr<'lua> {
    /// Create a new VCFExpr object. This object will read a VCF file, evaluate a set of expressions.
    /// The expressions should return a boolean. Evaluations will stop on the first true expression.
    /// If a template is provided, the template will be evaluated in the same scope as the expression and used
    /// to generate the text output. If no template is provided, the VCF record will be written to the output.
    pub fn new(
        lua: &'lua Lua,
        vcf_path: String,
        expression: Vec<String>,
        template: Option<String>,
        lua_code: Vec<String>,
        lua_prelude: Option<String>,
        output: Option<String>,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        for path in lua_code {
            let code = std::fs::read_to_string(&path)?;
            match lua.load(&code).set_name(&path).exec() {
                Ok(_) => (),
                Err(e) => {
                    log::error!("Error loading Lua code from {}: {}", path, e);
                    return Err(e.into());
                }
            }
        }

        lua.load(crate::pprint::PPRINT).set_name("pprint").exec()?;

        let reader = match vcf_path.as_str() {
            "-" | "stdin" => bcf::Reader::from_stdin()?,
            _ => bcf::Reader::from_path(&vcf_path)?,
        };
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

        let header = bcf::header::Header::from_template(reader.header());

        let writer = if template.is_none() {
            EitherWriter::Vcf(if let Some(output) = output {
                let format = get_vcf_format(&output);
                let mut wtr =
                    bcf::Writer::from_path(&output, &header, !output.ends_with(".gz"), format)?;
                _ = wtr.set_threads(2);
                //let header_t = unsafe { rust_htslib::htslib::bcf_hdr_dup(reader.header().inner) };
                //hv = bcf::header::HeaderView::new(header_t);
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
            vcf_reader: reader,
            template,
            writer,
            expressions: exps,
            globals,
            variants_evaluated: 0,
            variants_passing: 0,
        })
    }

    /// Return a reference to the bcf::Reader object.
    pub fn reader(&mut self) -> &mut bcf::Reader {
        &mut self.vcf_reader
    }

    pub fn evaluate(&mut self, record: bcf::Record) -> std::io::Result<StringOrBool> {
        let mut variant = Variant::new(record);
        self.variants_evaluated += 1;
        let eval_result = self.lua.scope(|scope| {
            let ud = match scope.create_any_userdata_ref_mut(&mut variant) {
                Ok(ud) => ud,
                Err(e) => return Err(e),
            };
            match self.globals.raw_set("variant", ud) {
                Ok(_) => (),
                Err(e) => return Err(e),
            }
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
                                Ok(res) => Ok(StringOrBool::String(res)),
                                Err(e) => {
                                    log::error!("Error in template: {}", e);
                                    return Err(e);
                                }
                            };
                        }
                        return Ok(StringOrBool::Bool(true));
                    }
                    Ok(false) => {}
                }
            }

            Ok(StringOrBool::Bool(false))
        });
        match eval_result {
            Ok(b) => Ok(b),
            Err(e) => Err(std::io::Error::new(std::io::ErrorKind::Other, e)),
        }
    }

    pub fn write(&mut self, record: &bcf::Record, sob: &StringOrBool) -> std::io::Result<()> {
        match sob {
            StringOrBool::Bool(false) => Ok(()),
            StringOrBool::Bool(true) => {
                if let EitherWriter::Vcf(wtr) = &mut self.writer {
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
            StringOrBool::String(s) => match &mut self.writer {
                EitherWriter::Vcf(_wtr) => Err(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    "did not VCF writer with template",
                )),
                EitherWriter::File(f) => writeln!(f, "{}", s),
                EitherWriter::Stdout(f) => writeln!(f, "{}", s),
            },
        }
    }
}
