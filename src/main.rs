use clap::{Parser, Subcommand};

use std::io::Write;

use mlua::Lua;
use rust_htslib::bcf::{self, Read};

use vcfexpr::register;
use vcfexpr::variant::Variant;

use log::info;

/// Args take the arguments for clap.
/// Accept the path to VCF or BCF and the lua expressions
#[derive(Parser)]
#[command(version, about, author)]
struct Cli {
    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand)]
pub enum Commands {
    /// Filter the VCF or BCF file and optionally apply a template.
    /// If no template is given the output will be VCF/BCF
    Filter {
        /// Path to input VCF or BCF file
        path: String,

        /// boolean Lua expression(s) to filter the VCF or BCF file
        #[arg(short, long)]
        expression: Vec<String>,

        /// template expression in luau: https://luau-lang.org/syntax#string-interpolation. e.g. '{variant.chrom}:{variant.pos}'
        #[arg(short, long)]
        template: Option<String>,

        /// File(s) containing lua code to run. Can contain functions that will be called by the expressions.
        #[arg(short, long)]
        lua: Vec<String>,

        /// File containing lua code to run once before any variants are processed.
        #[arg(short = 'p', long)]
        lua_prelude: Option<String>,

        /// optional output file. Default is stdout.
        #[arg(short, long)]
        output: Option<String>,
    },
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

fn get_vcf_format(path: &str) -> bcf::Format {
    if path.ends_with(".bcf") || path.ends_with(".bcf.gz") {
        bcf::Format::Bcf
    } else {
        bcf::Format::Vcf
    }
}

enum EitherWriter {
    Vcf(bcf::Writer),
    File(std::io::BufWriter<std::fs::File>),
    Stdout(std::io::BufWriter<std::io::Stdout>),
}
fn filter_main(
    path: String,
    expression: Vec<String>,
    template: Option<String>,
    lua_code: Vec<String>,
    lua_prelude: Option<String>,
    output: Option<String>,
) -> Result<(), Box<dyn std::error::Error>> {
    env_logger::init();
    let lua = Lua::new();

    for path in lua_code {
        let code = std::fs::read_to_string(&path)?;
        match lua.load(&code).set_name(path).exec() {
            Ok(_) => {}
            Err(e) => {
                log::error!("error in lua code: {}", e);
                return Err(e.into());
            }
        }
    }
    lua.load(vcfexpr::pprint::PPRINT)
        .set_name("pprint")
        .exec()?;

    //  open the VCF or BCF file
    let mut reader = bcf::Reader::from_path(path)?;
    _ = reader.set_threads(2);
    // create a new header from the reader
    //let header_t = unsafe { rust_htslib::htslib::bcf_hdr_dup(reader.header().inner) };

    //let mut hv = bcf::header::HeaderView::new(header_t);

    register(&lua)?;
    let globals = lua.globals();
    let template = process_template(template, &lua);

    let exps: Vec<_> = expression
        .iter()
        .map(|exp| {
            lua.load(exp)
                .set_name(exp)
                .into_function()
                .expect("error in expression")
        })
        .collect();

    let mut passing = 0;
    let mut hv =
        bcf::header::HeaderView::new(unsafe { rust_htslib::htslib::bcf_hdr_dup(reader.header().inner)});

    if let Some(lua_code) = lua_prelude {
        let code = std::fs::read_to_string(lua_code)?;
        lua.scope(|scope| {
            globals.raw_set("header", scope.create_any_userdata_ref_mut(&mut hv)?)?;
            lua.load(&code).exec()
        })?;
    }
    let header = bcf::header::Header::from_template(reader.header());

    let mut writer = if template.is_none() {
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
    let header_t = unsafe { rust_htslib::htslib::bcf_hdr_dup(reader.header().inner) };
    let hv = bcf::header::HeaderView::new(header_t);

    loop {
        let mut variant = reader.empty_record();
        match reader.read(&mut variant) {
            Some(Ok(_)) => {} //unsafe {
                //_ = rust_htslib::htslib::bcf_subset_format(hv.inner, variant.inner)
            //},
            Some(Err(e)) => return Err(e.into()),
            None => break,
        }
        if let EitherWriter::Vcf(ref mut w) = writer {
            w.translate(&mut variant);
        }
        let mut variant = Variant::new(variant);
        match check_variant(
            &lua,
            &mut variant,
            &hv,
            &exps,
            &template,
            &globals,
            &mut writer,
        ) {
            Ok(true) => {
                passing += 1;
            }
            Ok(false) => {}
            Err(e) => return Err(e.into()),
        }
    }
    info!("passing variants: {}", passing);
    Ok(())
}

enum StringOrBool {
    String(String),
    Bool(bool),
}

fn check_variant(
    lua: &Lua,
    variant: &mut Variant,
    header: &bcf::header::HeaderView,
    exps: &Vec<mlua::Function<'_>>,
    template: &Option<mlua::Function<'_>>,
    globals: &mlua::Table<'_>,
    writer: &mut EitherWriter,
) -> mlua::Result<bool> {
    let mut any_expression_true: bool = false;
    let result = lua.scope(|scope| {
        // TODO: ref_mut to allow setting.
        globals.raw_set("header", scope.create_any_userdata_ref(header)?)?;
        globals.raw_set("variant", scope.create_any_userdata_ref_mut(variant)?)?;

        for exp in exps {
            let result = exp.call::<_, bool>(());
            match result {
                Err(e) => return Err(e),
                Ok(true) => {
                    any_expression_true = true;
                    if let Some(template) = template {
                        return match template.call::<_, String>(()) {
                            Ok(s) => Ok(StringOrBool::String(s)),
                            Err(e) => {
                                log::error!("error in template: {}", e);
                                Err(e)
                            }
                        };
                    } else {
                        return Ok(StringOrBool::Bool(true));
                    }
                }
                Ok(false) => {}
            }
        }
        Ok(StringOrBool::Bool(false))
    });

    match result {
        Err(e) => Err(e),
        Ok(StringOrBool::Bool(false)) => Ok(false),
        Ok(StringOrBool::Bool(true)) => {
            match writer {
                EitherWriter::Vcf(w) => {
                    w.write(variant.record()).expect("error writing variant");
                }
                EitherWriter::File(_) | EitherWriter::Stdout(_) => {
                    return Err(mlua::Error::RuntimeError(
                        "File/Stdouput output not supported for variant".to_string(),
                    ));
                }
            }
            Ok(true)
        }
        Ok(StringOrBool::String(s)) => {
            match writer {
                EitherWriter::Vcf(_w) => {
                    // return a new error with msg
                    let msg = "VCF output not supported for template";
                    return Err(mlua::Error::RuntimeError(msg.to_string()));
                }
                EitherWriter::File(w) => {
                    writeln!(w, "{}", s).expect("error writing variant");
                }
                EitherWriter::Stdout(w) => {
                    writeln!(w, "{}", s).expect("error writing variant");
                }
            };
            Ok(true)
        }
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Cli::parse();
    match args.command {
        Some(Commands::Filter {
            path,
            expression,
            template,
            lua: lua_code,
            lua_prelude,
            output,
        }) => {
            filter_main(path, expression, template, lua_code, lua_prelude, output)?;
        }
        None => {
            println!("No command provided");
        }
    }
    Ok(())
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
