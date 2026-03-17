//! Refactored VCFExpress using the ScriptEngine trait
//!
//! This is a demonstration of how VCFExpress would look after being refactored
//! to use the ScriptEngine abstraction layer.

use rust_htslib::bcf::{
    self,
    header::{TagLength, TagType},
    Read,
};
use std::{collections::HashMap, hash::Hash, io::Write};

use crate::variant::{HeaderMap, Variant};
use crate::script_engine::{ScriptEngine, CompiledExpression, CompiledTemplate};

/// Refactored VCFExpress that uses the ScriptEngine trait
pub struct VCFExpress {
    /// The scripting engine (Lua, JavaScript, etc.)
    engine: Box<dyn ScriptEngine>,
    /// VCF reader
    vcf_reader: Option<bcf::Reader>,
    /// Compiled template for output formatting
    template: Option<CompiledTemplate>,
    /// Output writer
    writer: Option<EitherWriter>,
    /// Compiled filter expressions
    expressions: Vec<CompiledExpression>,
    /// Compiled expressions for setting INFO fields
    set_expressions: HashMap<InfoFormat, ((TagType, TagLength), CompiledExpression)>,
    /// Count of evaluated variants
    variants_evaluated: usize,
    /// Count of variants that passed filters
    variants_passing: usize,
}

/// `StringOrVariant` allows `evaluate` to return either a string, an owned VCF record, or nothing.
pub enum StringOrVariant {
    String(String),
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
    Format(String),
}

impl VCFExpress {
    /// Create a new VCFExpress object with the specified scripting engine
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

        // Compile set expressions
        let mut set_expressions = HashMap::new();
        // TODO: Implement set_expressions parsing and compilation

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
        } else if let Some(output) = output {
            EitherWriter::File(std::io::BufWriter::new(std::fs::File::create(output)?))
        } else {
            EitherWriter::Stdout(std::io::BufWriter::new(std::io::stdout()))
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

    /// Get a reference to the VCF reader
    pub fn reader(&self) -> &bcf::Reader {
        self.vcf_reader.as_ref().unwrap()
    }

    /// Get a mutable reference to the writer
    pub fn writer(&mut self) -> &mut EitherWriter {
        self.writer.as_mut().unwrap()
    }

    /// Evaluate a variant against all expressions
    pub fn evaluate(
        &mut self,
        mut record: bcf::Record,
        header: &HeaderMap,
        _header_map: HeaderMap,
    ) -> Result<StringOrVariant, Box<dyn std::error::Error>> {
        self.variants_evaluated += 1;

        // Create a Variant object from the record
        let variant = Variant::new(&record);

        // Evaluate all filter expressions
        let mut passes = false;
        for expr in &self.expressions {
            if self.engine.evaluate_variant(&variant, expr)? {
                passes = true;
                break; // Stop at first true expression
            }
        }

        if !passes {
            return Ok(StringOrVariant::None);
        }

        self.variants_passing += 1;

        // Apply set expressions if any
        for (field_info, expr) in &self.set_expressions {
            let field = match field_info {
                InfoFormat::Info(f) => f,
                InfoFormat::Format(_) => continue, // TODO: Handle FORMAT fields
            };
            self.engine.set_info_field(&mut record, field, expr)?;
        }

        // If there's a template, render it
        if let Some(template) = &self.template {
            let rendered = self.engine.render_template(&variant, template)?;
            Ok(StringOrVariant::String(rendered))
        } else {
            // Return the variant as-is for VCF output
            Ok(StringOrVariant::Variant(Some(record)))
        }
    }

    /// Get statistics
    pub fn stats(&self) -> (usize, usize) {
        (self.variants_passing, self.variants_evaluated)
    }
}