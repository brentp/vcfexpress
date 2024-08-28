use rusty_v8 as v8;
use rust_htslib::bcf::{
    self,
    header::{TagLength, TagType},
    Read,
};
use std::{collections::HashMap, hash::Hash, io::Write};

use crate::variant::{HeaderMap, Variant};

/// VCFExpress is the only entry-point for this library.
pub struct VCFExpress {
    isolate: v8::OwnedIsolate,
    context: v8::Global<v8::Context>,
    vcf_reader: Option<bcf::Reader>,
    template: Option<v8::Global<v8::Function>>,
    writer: Option<EitherWriter>,
    expressions: Vec<v8::Global<v8::Function>>,
    set_expressions: HashMap<InfoFormat, ((TagType, TagLength), v8::Global<v8::Function>)>,
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

/*
fn process_template(template: Option<String>, isolate: &v8::OwnedIsolate) -> Option<v8::Global<v8::Function>> {
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
        let mut handle_scope = v8::HandleScope::new(isolate);
        let context = v8::Context::new(&mut handle_scope);
        let global = context.global(&mut handle_scope);
        let source = v8::String::new(&mut handle_scope, &expr).unwrap();
        let script = v8::Script::compile(&mut handle_scope, source, None).unwrap();
        let function = script.run(&mut handle_scope).unwrap().into();
        Some(v8::Global::new(&mut handle_scope, function))
    } else {
        None
    }
}
    */

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

impl VCFExpress {
    /// Create a new VCFExpress object. This object will read a VCF file, evaluate a set of expressions.
    /// The expressions should return a boolean. Evaluations will stop on the first true expression.
    /// If a template is provided, the template will be evaluated in the same scope as the expression and used
    /// to generate the text output. If no template is provided, the VCF record will be written to the output.
    /// The template is a [luau string template].
    ///
    /// [luau string template]: https://luau.org/syntax#string-interpolation
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        vcf_path: String,
        expression: Vec<String>,
        set_expression: Vec<String>,
        template: Option<String>,
        js_prelude: Vec<String>,
        output: Option<String>,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        // Initialize V8
        let platform = v8::new_default_platform(0, false).make_shared();
        v8::V8::initialize_platform(platform);
        v8::V8::initialize();

        let mut isolate = v8::Isolate::new(Default::default());
        let context = {
            let handle_scope = &mut v8::HandleScope::new(&mut isolate);
            v8::Context::new(handle_scope)
        };
        let scope = &mut v8::HandleScope::with_context(&mut isolate, &context);
        let global = context.global(scope);

        // Register VCF functions and objects
        // This part needs to be implemented to expose VCF functionality to JavaScript

        // Compile expressions
        let expressions = expression.iter().map(|exp| {
            let source = v8::String::new(scope, exp).unwrap();
            let script = v8::Script::compile(scope, source, None).unwrap();
            v8::Global::new(scope, script.run(scope).unwrap().into())
        }).collect();

        // Similar changes for set_expressions and template

        // ... rest of the implementation

        let vcf_reader = bcf::Reader::from_path(&vcf_path)?;
        let header = vcf_reader.header().clone();

        Ok(VCFExpress {
            isolate,
            context: v8::Global::new(scope, context),
            vcf_reader: Some(vcf_reader),
            template: None, //process_template(template, &isolate),
            writer: Some(EitherWriter::Vcf(bcf::Writer::from_path(&output.unwrap_or_else(|| "-".to_string()), &header, true, bcf::Format::Vcf)?)),
            expressions,
            set_expressions: HashMap::new(), // Initialize this properly
            variants_evaluated: 0,
            variants_passing: 0,
        })
    }

    /// Run the code in the luau sandboxed environment.
    /// https://luau.org/sandbox
    pub fn sandbox(&mut self, _sandbox: bool) -> Result<(), Box<dyn std::error::Error>> {
        // Implement sandbox logic for V8 here
        // This is a placeholder, you'll need to implement actual V8 sandboxing
        Ok(())
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

    // this is called from in the scope and lets us evaluate the info expressions.
    // we collect the results to be used outside the scope where we can get a mutable variant.

    pub fn evaluate_info_expressions(
        &mut self,
        info_results: &mut HashMap<String, InfoFormatValue>,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut scope = v8::HandleScope::with_context(&mut self.isolate, &self.context);
        
        for (inf, ((tagtyp, _taglen), expr)) in self.set_expressions.iter() {
            if let InfoFormat::Info(tag) = inf {
                let function = v8::Local::new(&mut scope, expr);
                let global = scope.get_current_context().global(&mut scope);
                let result = function.call(&mut scope, global.into(), &[]);
                
                if let Some(result) = result {
                    let t = match tagtyp {
                        TagType::Flag => {
                            let b = result.boolean_value(&mut scope);
                            InfoFormatValue::Bool(b)
                        }
                        TagType::Float => {
                            let f = result.number_value(&mut scope).unwrap_or(0.0) as f32;
                            InfoFormatValue::Float(f)
                        }
                        TagType::Integer => {
                            let i = result.integer_value(&mut scope).unwrap_or(0) as i32;
                            InfoFormatValue::Integer(i)
                        }
                        TagType::String => {
                            let s = result.to_string(&mut scope).unwrap().to_rust_string_lossy(&mut scope);
                            InfoFormatValue::String(s)
                        }
                    };
                    info_results.insert(tag.clone(), t);
                }
            }
        }
        Ok(())
    }

    /// Evaluate the expressions and optional template for a single record.
    pub fn evaluate(
        &mut self,
        record: bcf::Record,
        header_map: HeaderMap,
    ) -> std::io::Result<StringOrVariant> {
        let mut variant = Variant::new(record, header_map);
        self.variants_evaluated += 1;

        let mut scope = v8::HandleScope::with_context(&mut self.isolate, &self.context);
        let global = scope.get_current_context().global(&mut scope);

        // Create JavaScript Variant object
        let variant_obj = v8::ObjectTemplate::new(&mut scope);
        variant_obj.set_internal_field_count(1);
        let variant_instance = variant_obj.new_instance(&mut scope).unwrap();
        variant_instance.set_internal_field(0, v8::External::new(&mut scope, &variant as *const _ as *mut std::ffi::c_void).into());

        global.set(
            &mut scope,
            v8::String::new(&mut scope, "variant").unwrap().into(),
            variant_instance.into(),
        ).unwrap();

        let mut result = StringOrVariant::None;

        for exp in &self.expressions {
            let function = v8::Local::new(&mut scope, exp);
            let undefined = v8::undefined(&mut scope);
            let global_context = v8::Local::new(&mut scope, self.context);
            let result_value = function.call(&mut scope, global_context.into(), &[]);

            if let Some(result_value) = result_value {
                if result_value.is_true() {
                    self.variants_passing += 1;
                    if let Some(template) = &self.template {
                        let template_function = v8::Local::new(&mut scope, template);
                        let template_result = template_function.call(&mut scope, undefined.into(), &[]);
                        if let Some(template_result) = template_result {
                            if template_result.is_string() {
                                let string = template_result.to_string(&mut scope).unwrap();
                                result = StringOrVariant::String(string.to_rust_string_lossy(&mut scope));
                            }
                        }
                    } else {
                        result = StringOrVariant::Variant(Some(variant.take()));
                    }
                    break;
                }
            }
        }

        // ... rest of the implementation

        Ok(result)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use rusty_v8 as v8;

    fn setup_v8() -> v8::OwnedIsolate {
        let platform = v8::new_default_platform(0, false).make_shared();
        v8::V8::initialize_platform(platform);
        v8::V8::initialize();
        v8::Isolate::new(Default::default())
    }

    /*
    #[test]
    fn test_process_template_with_none() {
        let isolate = setup_v8();
        assert_eq!(process_template(None, &isolate), None);
    }

    #[test]
    fn test_process_template_with_backticks() {
        let isolate = setup_v8();
        let template = Some("`console.log('Hello, World!')`".to_string());
        let result = process_template(template, &isolate);
        assert!(result.is_some());
    }

    #[test]
    fn test_process_template_without_backticks() {
        let isolate = setup_v8();
        let template = Some("console.log('Hello, World!')".to_string());
        let result = process_template(template, &isolate);
        assert!(result.is_some());
    }

    #[test]
    fn test_process_template_with_return() {
        let mut isolate = setup_v8();
        let template = Some("return `42`".to_string());
        let result = process_template(template, &isolate);
        assert!(result.is_some());

        let scope = &mut v8::HandleScope::new(&mut isolate);
        let context = v8::Context::new(scope);
        let scope = &mut v8::ContextScope::new(scope, context);

        let function = v8::Local::new(scope, result.unwrap());
        let result = function.call(scope, context.global(scope).into(), &[]);
        
        assert!(result.is_some());
        let result = result.unwrap();
        assert!(result.is_string());
        let result_str = result.to_string(scope).unwrap().to_rust_string_lossy(scope);
        assert_eq!(result_str, "42");
    }

    #[test]
    #[should_panic(expected = "SyntaxError")]
    fn test_process_template_with_invalid_js() {
        let isolate = setup_v8();
        let template = Some("return [invalid_js_code".to_string());
        process_template(template, &isolate);
    }
    */
}