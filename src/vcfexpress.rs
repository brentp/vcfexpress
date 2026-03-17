use rust_htslib::bcf::{
    self,
    header::{TagLength, TagType},
    Read,
};
use std::{collections::HashMap, hash::Hash, io::Write};

use crate::variant::{HeaderMap, Variant};
use crate::script_engine::{ScriptEngine, CompiledExpression, CompiledTemplate};

#[cfg(feature = "lua")]
use mlua::Lua;

/// VCFExpress is the only entry-point for this library.
pub struct VCFExpress {
    /// The scripting engine (Lua, JavaScript, etc.)
    engine: Box<dyn ScriptEngine>,
    vcf_reader: Option<bcf::Reader>,
    template: Option<CompiledTemplate>,
    writer: Option<EitherWriter>,
    expressions: Vec<CompiledExpression>,
    set_expressions: HashMap<InfoFormat, ((TagType, TagLength), CompiledExpression)>,
    variants_evaluated: usize,
    variants_passing: usize,
}

/// `StringOrVariant` allows `evaluate` to return either a string, an owned VCF record, or nothing.
pub enum StringOrVariant {
    String(String),
    // Variant(None) is used since we sometimes can't take ownership of the
    // bcf::Record right away so we set Variant(None) and later replace
    // with Variant(Some(Record)).
    Variant(Option<bcf::Record>),
    None,
}

/// `EitherWriter` encapsulates the different types of writers we can use.
/// `File` and `Stdout` are for template output and `Vcf` is for VCF records.
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

    pub fn write(&mut self, sob: &mut StringOrVariant) -> std::io::Result<u32> {
        match sob {
            StringOrVariant::None => Ok(0),
            StringOrVariant::Variant(None) => Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                "expected VCF record got None",
            )),
            StringOrVariant::Variant(Some(ref mut record)) => {
                if let EitherWriter::Vcf(ref mut wtr) = self {
                    match wtr.write(record) {
                        Ok(_) => Ok(1),
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
                EitherWriter::File(ref mut f) => writeln!(f, "{}", s).map(|_| 1),
                EitherWriter::Stdout(ref mut f) => writeln!(f, "{}", s).map(|_| 1),
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

#[derive(Debug, PartialEq, Eq, Hash)]
enum InfoFormat {
    Info(String),
    #[allow(dead_code)]
    Format(String),
}

impl VCFExpress {
    /// Create a new VCFExpress object using a ScriptEngine.
    /// This object will read a VCF file, evaluate a set of expressions.
    /// The expressions should return a boolean. Evaluations will stop on the first true expression.
    /// If a template is provided, the template will be evaluated in the same scope as the expression and used
    /// to generate the text output. If no template is provided, the VCF record will be written to the output.
    /// The template syntax depends on the scripting engine (Luau for Lua, template literals for JavaScript).
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_engine(
        engine: Box<dyn ScriptEngine>,
        vcf_path: String,
        expression: Vec<String>,
        set_expression: Vec<String>,
        template: Option<String>,
        prelude_files: Vec<String>,
        output: Option<String>,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let mut reader = match vcf_path.as_str() {
            "-" | "stdin" => bcf::Reader::from_stdin()?,
            _ => bcf::Reader::from_path(&vcf_path)?,
        };
        _ = reader.set_threads(2);

        // Create a mutable reference to the engine for configuration
        let mut engine_mut = engine;

        // Register types with the scripting engine
        let header_view = bcf::header::HeaderView::new(unsafe {
            rust_htslib::htslib::bcf_hdr_dup(reader.header().inner)
        });
        engine_mut.register_types(&header_view)?;

        // Load all prelude files
        for prelude_path in prelude_files {
            let prelude_code = std::fs::read_to_string(&prelude_path)?;
            engine_mut.load_prelude(&prelude_code)?;
        }

        // Compile all filter expressions
        let mut expressions = Vec::new();
        for expr in expression {
            let compiled = engine_mut.compile_expression(&expr)?;
            expressions.push(compiled);
        }

        // Compile the template if provided
        let template = if let Some(tpl) = template {
            Some(engine_mut.compile_template(&tpl)?)
        } else {
            None
        };

        // Parse and compile set expressions
        let mut set_expressions = HashMap::new();
        for exp in set_expression {
            let name_exp = exp
                .split_once('=')
                .expect("invalid info expression should have name=$expression");
            let field_name = name_exp.0.to_string();
            let expression_str = name_exp.1;

            // Get the field type from header
            let tag_type = header_view
                .info_type(field_name.as_bytes())
                .unwrap_or_else(|_| {
                    panic!("ERROR: info field '{}' not found. Make sure it was added to the header in prelude if needed.", field_name)
                });

            let compiled_expr = engine_mut.compile_expression(expression_str)?;
            set_expressions.insert(
                InfoFormat::Info(field_name),
                ((tag_type.0, tag_type.1), compiled_expr),
            );
        }

        let header = bcf::header::Header::from_template(&header_view);

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

        Ok(VCFExpress {
            engine: engine_mut,
            vcf_reader: Some(reader),
            template,
            writer: Some(writer),
            expressions,
            set_expressions,
            variants_evaluated: 0,
            variants_passing: 0,
        })
    }

    /// Legacy constructor using Lua for backward compatibility
    #[cfg(feature = "lua")]
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        lua: Lua,
        vcf_path: String,
        expression: Vec<String>,
        set_expression: Vec<String>,
        template: Option<String>,
        lua_prelude: Vec<String>,
        output: Option<String>,
        sandbox: bool,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        // Create a temporary VCFExpress to handle the legacy initialization
        // This is a hack to maintain backward compatibility
        // In the future, users should migrate to new_with_engine
        panic!("Legacy constructor not yet implemented. Use new_with_engine instead.");
    }

    /// Note: The sandbox method is now handled by the ScriptEngine configuration
    /// This method is kept for compatibility but does nothing
    pub fn sandbox(&mut self, _sandbox: bool) -> Result<(), Box<dyn std::error::Error>> {
        // Sandbox is now configured when creating the ScriptEngine
        Ok(())
    }

    #[allow(clippy::type_complexity)]
    fn load_info_expressions(
        lua: &Lua,
        hv: &mut bcf::header::HeaderView,
        info_expressions: Vec<String>,
    ) -> Result<
        HashMap<InfoFormat, ((TagType, TagLength), mlua::Function)>,
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

    /// Add code to the script engine. This code will be available to the expressions and the template.
    /// These are not the variant expressions, but rather additional code that can be used as a library.
    pub fn add_code(&mut self, path: &str) -> Result<(), Box<dyn std::error::Error>> {
        let code = std::fs::read_to_string(path)
            .map_err(|e| format!("Error reading file {}: {}", path, e))?;
        self.engine.load_prelude(&code)?;
        Ok(())
    }

    /// Add lua code to the Lua interpreter. This code will be available to the expressions and the template.
    /// These are not the variant expressions, but rather additional Lua code that can be used as a library.
    /// Deprecated: Use add_code instead.
    #[deprecated(note = "Use add_code instead")]
    pub fn add_lua_code(&mut self, path: &str) -> Result<(), Box<dyn std::error::Error>> {
        self.add_code(path)
    }

    /// Take ownership of the the bcf::Reader object.
    /// This must be called before using `evaluate`
    pub fn reader(&mut self) -> bcf::Reader {
        self.vcf_reader.take().expect("reader already taken")
    }

    /// Take ownership of the the Writer enum.
    /// This must be called before using `evaluate`
    pub fn writer(&mut self) -> EitherWriter {
        self.writer.take().expect("writer already taken")
    }

    /// Evaluate the expressions and optional template for a single record using the ScriptEngine.
    pub fn evaluate(
        &mut self,
        record: bcf::Record,
        header: &bcf::header::HeaderView,
        header_map: HeaderMap,
    ) -> std::io::Result<StringOrVariant> {
        let mut variant = Variant::new(record, header_map.clone());
        self.variants_evaluated += 1;

        // Evaluate all filter expressions - stop at first true
        let mut passes = false;
        for expr in &self.expressions {
            match self.engine.evaluate_variant(&variant, expr) {
                Ok(true) => {
                    passes = true;
                    break;
                }
                Ok(false) => continue,
                Err(e) => {
                    log::error!("Error evaluating expression: {}", e);
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::Other,
                        format!("Script evaluation error: {}", e),
                    ));
                }
            }
        }

        if !passes {
            return Ok(StringOrVariant::None);
        }

        self.variants_passing += 1;

        // Extract the record before potentially applying set expressions
        let mut record = variant.take();

        // Apply set expressions if any
        for (field_info, ((_, _), expr)) in &self.set_expressions {
            let field = match field_info {
                InfoFormat::Info(f) => f,
                InfoFormat::Format(_) => continue, // TODO: Handle FORMAT fields
            };

            if let Err(e) = self.engine.set_info_field(&mut record, field, expr) {
                log::error!("Error setting info field '{}': {}", field, e);
                return Err(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!("Failed to set info field: {}", e),
                ));
            }
        }

        // If there's a template, render it
        if let Some(template) = &self.template {
            // Need to recreate variant for template rendering
            let variant_for_template = Variant::new(record, header_map);
            match self.engine.render_template(&variant_for_template, template) {
                Ok(rendered) => Ok(StringOrVariant::String(rendered)),
                Err(e) => {
                    log::error!("Error rendering template: {}", e);
                    Err(std::io::Error::new(
                        std::io::ErrorKind::Other,
                        format!("Template rendering error: {}", e),
                    ))
                }
            }
        } else {
            // Return the variant as-is for VCF output
            Ok(StringOrVariant::Variant(Some(record)))
        }
    }

    /// Legacy evaluate method for backward compatibility
    /// Now delegates to the new ScriptEngine-based evaluate method
    #[deprecated(note = "Use evaluate method instead")]
    pub fn evaluate_legacy(
        &mut self,
        record: bcf::Record,
        header: &bcf::header::HeaderView,
        header_map: HeaderMap,
    ) -> std::io::Result<StringOrVariant> {
        self.evaluate(record, header, header_map)
    }
}

