use clap::{Parser, Subcommand};
use rust_htslib::bcf::record::{Genotype, GenotypeAllele};
use std::borrow::{Borrow, BorrowMut};
use std::io::Write;

use env_logger;
use mlua::prelude::LuaValue;
use mlua::{AnyUserData, Lua, MetaMethod, UserDataFields, UserDataMethods, Value};
use rust_htslib::bcf::{self, Read};
use std::sync::Arc;

use log::{debug, info, log_enabled, Level};

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
    /// If no template is given the output will be VCF/BCF (TODO)
    Filter {
        /// Path to input VCF or BCF file
        path: String,

        /// boolean Lua expression(s) to filter the VCF or BCF file
        #[arg(short, long)]
        expression: Vec<String>,

        /// template expression in luau: https://luau-lang.org/syntax#string-interpolation. e.g. '{variant.chrom}:{variant.pos}'
        #[arg(short, long)]
        template: Option<String>,

        /// optional output file. Default is stdout.
        #[arg(short, long)]
        output: Option<String>,
    },
}

struct Variant(bcf::Record);

fn register_variant(lua: &Lua) -> mlua::Result<()> {
    lua.register_userdata_type::<Variant>(|reg| {
        reg.add_meta_function(
            MetaMethod::Index,
            |_lua, (_, name): (AnyUserData, String)| {
                let msg = format!("field '{}' variant.{} not found", name, name);
                Err::<LuaValue<'_>, mlua::Error>(mlua::Error::RuntimeError(msg))
            },
        );
        reg.add_field_method_get("chrom", |_, this: &Variant| {
            let c = this
                .0
                .rid()
                .map(|id| this.0.header().rid2name(id))
                .unwrap_or(Ok(b""))
                .map(|c| unsafe { String::from_utf8_unchecked(c.to_vec()) });
            c.map_err(|e| mlua::Error::ExternalError(Arc::new(e)))
        });
        reg.add_field_method_get("qual", |_, this: &Variant| Ok(this.0.qual()));
        reg.add_field_method_set("qual", |_, this: &mut Variant, val: f32| {
            this.0.set_qual(val);
            Ok(())
        });

        reg.add_field_method_get("pos", |_, this: &Variant| Ok(this.0.pos()));
        reg.add_field_method_set("pos", |_, this: &mut Variant, val: i64| {
            this.0.set_pos(val as i64);
            Ok(())
        });
        reg.add_field_method_get("filters", |lua: &Lua, this: &Variant| {
            let f = this.0.filters();
            let t = lua.create_table().expect("error creating table");
            let h = this.0.header();
            for (i, id) in f.into_iter().enumerate() {
                let filter = unsafe { String::from_utf8_unchecked(h.id_to_name(id)) };
                t.raw_set(i + 1, filter).expect("error setting value");
            }
            Ok(Value::Table(t))
        });
        reg.add_field_method_set(
            "filters",
            |_, this: &mut Variant, filter: String| match this
                .0
                .set_filters(&vec![filter.as_bytes()])
            {
                Err(e) => Err(mlua::Error::ExternalError(Arc::new(e))),
                Ok(_) => Ok(()),
            },
        );
        reg.add_field_method_get("id", |lua: &Lua, this: &Variant| {
            let id = this.0.id();
            Ok(Value::String(unsafe {
                lua.create_string(String::from_utf8_unchecked(id.to_vec()))?
            }))
        });
        reg.add_field_method_set(
            "id",
            |_lua: &Lua, this: &mut Variant, val: String| match this.0.set_id(val.as_bytes()) {
                Err(e) => Err(mlua::Error::ExternalError(Arc::new(e))),
                Ok(_) => Ok(()),
            },
        );

        reg.add_field_method_get("REF", |lua: &Lua, this: &Variant| {
            let ref_allele = this.0.alleles()[0];
            Ok(Value::String(unsafe {
                lua.create_string(String::from_utf8_unchecked(ref_allele.to_vec()))?
            }))
        });
        reg.add_field_method_set("REF", |_lua: &Lua, this: &mut Variant, val: String| {
            let mut alleles = this.0.alleles();
            alleles[0] = val.as_bytes();
            Ok(())
        });

        reg.add_field_method_get("ALT", |lua: &Lua, this: &Variant| {
            let alt_alleles = this.0.alleles();
            let count = alt_alleles.len() - 1;
            let t = lua
                .create_table_with_capacity(count, 0)
                .expect("error creating table");
            for (i, allele) in alt_alleles.iter().skip(1).enumerate() {
                t.raw_set(i + 1, unsafe {
                    String::from_utf8_unchecked(allele.to_vec())
                })
                .expect("error setting value");
            }
            if t.is_empty() {
                t.raw_set(1, lua.create_string(b".")?)
                    .expect("error setting value");
            }
            Ok(Value::Table(t))
        });
        reg.add_field_method_get("FILTER", |lua: &Lua, this: &Variant| {
            let f = this.0.filters();
            let h = this.0.header();
            for id in f.into_iter() {
                let filter = unsafe { String::from_utf8_unchecked(h.id_to_name(id)) };
                return Ok(Value::String(lua.create_string(&filter)?));
            }
            return Ok(Value::Nil);
        });

        reg.add_method("format", |lua: &Lua, this: &Variant, format: String| {
            let fmt = this.0.format(format.as_bytes());
            let typ = this.0.header().format_type(format.as_bytes());
            let (typ, num) = match typ {
                Err(e) => return Err(mlua::Error::ExternalError(Arc::new(e))),
                Ok(typ) => typ,
            };
            let n_samples = this.0.sample_count() as usize;
            //eprintln!("{:?} {:?}", typ, num);
            let t = lua
                .create_table_with_capacity(n_samples, 0)
                .expect("error creating table");
            return match typ {
                bcf::header::TagType::Integer => fmt
                    .integer()
                    .map(|v| {
                        if matches!(num, bcf::header::TagLength::Fixed(1)) {
                            for (i, vals) in v.iter().enumerate() {
                                t.raw_set(i + 1, vals[0]).expect("error setting value");
                            }
                        } else {
                            for (i, vals) in v.iter().enumerate() {
                                let ti = lua
                                    .create_table_with_capacity(vals.len(), 0)
                                    .expect("error creating table");
                                for (j, val) in vals.iter().enumerate() {
                                    ti.raw_set(j + 1, *val).expect("error setting value");
                                }
                                t.raw_set(i + 1, ti).expect("error setting value");
                            }
                        }
                        Ok::<LuaValue<'_>, mlua::Error>(Value::Table(t))
                    })
                    .map_err(|e| mlua::Error::ExternalError(Arc::new(e))),
                bcf::header::TagType::Float => fmt
                    .float()
                    .map(|v| {
                        if matches!(num, bcf::header::TagLength::Fixed(1)) {
                            for (i, vals) in v.iter().enumerate() {
                                t.raw_set(i + 1, vals[0]).expect("error setting value");
                            }
                        } else {
                            for (i, vals) in v.iter().enumerate() {
                                let ti = lua
                                    .create_table_with_capacity(vals.len(), 0)
                                    .expect("error creating table");
                                for (j, val) in vals.iter().enumerate() {
                                    ti.raw_set(j + 1, *val).expect("error setting value");
                                }
                                t.raw_set(i + 1, ti).expect("error setting value");
                            }
                        }
                        Ok::<LuaValue<'_>, mlua::Error>(Value::Table(t))
                    })
                    .map_err(|e| mlua::Error::ExternalError(Arc::new(e))),

                bcf::header::TagType::String => Ok(fmt.string().map_or_else(
                    |_e| {
                        if log_enabled!(Level::Debug) {
                            debug!("format tag {} not found", format);
                        }
                        Ok(Value::Nil)
                    },
                    |v| {
                        for (i, vals) in v.iter().enumerate() {
                            t.raw_set(i + 1, unsafe { String::from_utf8_unchecked(vals.to_vec()) })
                                .expect("error setting value");
                        }
                        Ok::<LuaValue<'_>, mlua::Error>(Value::Table(t))
                    },
                )),

                _ => unimplemented!("format type {:?}", typ),
            };
        });

        reg.add_method(
            "info",
            |lua: &Lua, this: &Variant, (key, index): (String, Option<usize>)| {
                let mut info = this.0.info(key.as_bytes()); /* only need mut for .flag */
                let typ = this.0.header().info_type(key.as_bytes());
                let (typ, num) = match typ {
                    Err(e) => return Err(mlua::Error::ExternalError(Arc::new(e))),
                    Ok(typ) => typ,
                };
                return match typ {
                    bcf::header::TagType::Integer => info
                        .integer()
                        .map(|v| match v {
                            Some(v) => match (num, index) {
                                (bcf::header::TagLength::Fixed(1), None) => {
                                    Ok::<LuaValue<'_>, mlua::Error>(Value::Integer(v[0]))
                                }
                                (_, Some(i)) => {
                                    Ok::<LuaValue<'_>, mlua::Error>(Value::Integer(v[i as usize]))
                                }

                                _ => {
                                    let t = lua.create_table().expect("error creating table");
                                    for (i, val) in v.iter().enumerate() {
                                        t.raw_set(i + 1, *val).expect("error setting value");
                                    }
                                    Ok::<LuaValue<'_>, mlua::Error>(Value::Table(t))
                                }
                            },
                            None => Ok(Value::Nil),
                        })
                        .map_err(|e| mlua::Error::ExternalError(Arc::new(e))),
                    bcf::header::TagType::Float => info
                        .float()
                        .map(|v| match v {
                            Some(v) => match (num, index) {
                                (bcf::header::TagLength::Fixed(1), None) => {
                                    Ok::<LuaValue<'_>, mlua::Error>(Value::Number(v[0] as f64))
                                }
                                (_, Some(i)) => Ok::<LuaValue<'_>, mlua::Error>(Value::Number(
                                    v[i as usize] as f64,
                                )),
                                _ => {
                                    let t = lua.create_table().expect("error creating table");
                                    for (i, val) in v.iter().enumerate() {
                                        t.raw_set(i + 1, *val as f64).expect("error setting value");
                                    }
                                    Ok::<LuaValue<'_>, mlua::Error>(Value::Table(t))
                                }
                            },
                            None => Ok(Value::Nil),
                        })
                        .map_err(|e| mlua::Error::ExternalError(Arc::new(e))),
                    bcf::header::TagType::String => info
                        .string()
                        .map(|v| match v {
                            Some(v) => match (num, index) {
                                (bcf::header::TagLength::Fixed(1), None) => {
                                    Ok::<LuaValue<'_>, mlua::Error>(Value::String(
                                        lua.create_string(unsafe {
                                            String::from_utf8_unchecked(v[0].to_vec())
                                        })?,
                                    ))
                                }
                                (_, Some(i)) => Ok::<LuaValue<'_>, mlua::Error>(Value::String(
                                    lua.create_string(unsafe {
                                        String::from_utf8_unchecked(v[i as usize].to_vec())
                                    })?,
                                )),
                                _ => {
                                    let t = lua.create_table().expect("error creating table");
                                    for (i, s) in v.iter().enumerate() {
                                        t.raw_set(i + 1, unsafe {
                                            String::from_utf8_unchecked(s.to_vec())
                                        })
                                        .expect("error setting value");
                                    }
                                    Ok::<LuaValue<'_>, mlua::Error>(Value::Table(t))
                                }
                            },
                            None => Ok(Value::Nil),
                        })
                        .map_err(|e| mlua::Error::ExternalError(Arc::new(e))),
                    bcf::header::TagType::Flag => info
                        .flag()
                        .map(|v| Ok::<LuaValue<'_>, mlua::Error>(Value::Boolean(v)))
                        .map_err(|e| mlua::Error::ExternalError(Arc::new(e))),
                };
            },
        );
    })
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
    output: Option<String>,
) -> Result<(), Box<dyn std::error::Error>> {
    env_logger::init();
    let lua = Lua::new();

    //  open the VCF or BCF file
    let mut reader = bcf::Reader::from_path(&path)?;
    _ = reader.set_threads(2);
    // create a new header from the reader
    let header = bcf::Header::from_template(reader.header());

    let mut writer = if template.is_none() {
        EitherWriter::Vcf(if let Some(output) = output {
            let format = get_vcf_format(&output);
            let mut wtr =
                bcf::Writer::from_path(&output, &header, !output.ends_with(".gz"), format)?;
            _ = wtr.set_threads(2);
            wtr
        } else {
            bcf::Writer::from_stdout(&header, true, bcf::Format::Vcf)?
        })
    } else {
        if output.is_none() || output.as_ref().unwrap() == "-" {
            EitherWriter::Stdout(std::io::BufWriter::new(std::io::stdout()))
        } else {
            let file = std::fs::File::create(output.unwrap())?;
            EitherWriter::File(std::io::BufWriter::new(file))
        }
    };

    let globals = lua.globals();
    register_variant(&lua)?;
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

    for variant in reader.records() {
        let mut variant = Variant(variant?);
        match check_variant(&lua, &mut variant, &exps, &template, &globals, &mut writer) {
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

fn check_variant(
    lua: &Lua,
    variant: &mut Variant,
    exps: &Vec<mlua::Function<'_>>,
    template: &Option<mlua::Function<'_>>,
    globals: &mlua::Table<'_>,
    writer: &mut EitherWriter,
) -> mlua::Result<bool> {
    lua.scope(|scope| {
        let ud = scope.create_any_userdata_ref(variant)?;
        globals.raw_set("variant", ud)?;

        for exp in exps {
            let result = exp.call::<_, bool>(());
            match result {
                Ok(false) => {}
                Ok(true) => {
                    if let Some(template) = template {
                        let result = template.call::<_, String>(());
                        match result {
                            Ok(result) => match writer {
                                EitherWriter::Vcf(w) => {
                                    w.write(&variant.0).expect("error writing variant");
                                }
                                EitherWriter::File(w) => {
                                    writeln!(w, "{}", result).expect("error writing variant");
                                }
                                EitherWriter::Stdout(w) => {
                                    writeln!(w, "{}", result).expect("error writing variant");
                                }
                            },
                            Err(e) => return Err(e.into()),
                        }
                    } else {
                        match writer {
                            EitherWriter::Vcf(w) => {
                                w.write(&variant.0).expect("error writing variant");
                            }
                            _ => panic!("expected VCF output file."),
                        }
                    }
                    return Ok(true);
                }
                Err(e) => return Err(e.into()),
            }
        }
        Ok(false)
    })
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Cli::parse();
    match args.command {
        Some(Commands::Filter {
            path,
            expression,
            template,
            output,
        }) => {
            filter_main(path, expression, template, output)?;
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
