use log::{error, info, warn};
use mlua::prelude::LuaValue;
use mlua::{AnyUserData, Lua, MetaMethod, UserDataFields, UserDataMethods, Value};
use parking_lot::Mutex;
use rust_htslib::bcf::{self};
use std::sync::Arc;

pub struct Variant(bcf::Record);

impl Variant {
    pub fn new(record: bcf::Record) -> Self {
        Variant(record)
    }
    pub fn record(&self) -> &bcf::Record {
        &self.0
    }
    pub fn header(&self) -> &bcf::header::HeaderView {
        self.0.header()
    }
    pub fn take(self) -> bcf::Record {
        self.0
    }
}

use log::{debug, log_enabled, Level};

pub fn register_variant(lua: &Lua) -> mlua::Result<()> {
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

        reg.add_field_method_get("start", |_, this: &Variant| Ok(this.0.pos()));
        reg.add_field_method_get("stop", |_, this: &Variant| Ok(this.0.end()));
        reg.add_field_method_get("pos", |_, this: &Variant| Ok(this.0.pos()));
        reg.add_field_method_set("pos", |_, this: &mut Variant, val: i64| {
            this.0.set_pos(val);
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
            |_, this: &mut Variant, filter: String| match this.0.set_filters(&[filter.as_bytes()]) {
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
            let mut alleles = vec![val.as_bytes()];
            let alt_alleles = this
                .0
                .alleles()
                .iter()
                .skip(1)
                .map(|&a| a.to_owned())
                .collect::<Vec<_>>();
            alleles.extend(alt_alleles.iter().map(|a| &a[..]));

            match this.0.set_alleles(&alleles) {
                Ok(_) => Ok(()),
                Err(e) => Err(mlua::Error::ExternalError(Arc::new(e))),
            }
        });
        reg.add_field_method_set("ALT", |_lua: &Lua, this: &mut Variant, val: Vec<String>| {
            let ref_allele = this.0.alleles()[0].to_owned();
            let mut alleles = vec![&ref_allele[..]];
            alleles.extend(val.iter().map(|a| a.as_bytes()));

            match this.0.set_alleles(&alleles) {
                Ok(_) => Ok(()),
                Err(e) => Err(mlua::Error::ExternalError(Arc::new(e))),
            }
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
            if let Some(filter) = f.into_iter().next() {
                let filter = unsafe { String::from_utf8_unchecked(h.id_to_name(filter)) };
                return Ok(Value::String(lua.create_string(&filter)?));
            }
            Ok(Value::Nil)
        });
        reg.add_field_method_set(
            "FILTER",
            |_lua, this: &mut Variant, filter: String| match this
                .0
                .set_filters(&[filter.as_bytes()])
            {
                Err(e) => Err(mlua::Error::ExternalError(Arc::new(e))),
                Ok(_) => Ok(()),
            },
        );
        reg.add_field_method_get("genotypes", |_lua: &Lua, this: &Variant| {
            let genotypes = this.0.format(b"GT");
            match genotypes.integer() {
                Ok(genotypes) => {
                    let sb = crate::genotypes::Genotypes(Arc::new(Mutex::new(
                        crate::genotypes::I32Buffer(genotypes),
                    )));
                    Ok(sb)
                }
                Err(e) => Err(mlua::Error::RuntimeError(e.to_string())),
            }
        });

        reg.add_method("format", |lua: &Lua, this: &Variant, format: String| {
            let fmt = this.0.format(format.as_bytes());
            let typ = this.0.header().format_type(format.as_bytes());
            let (typ, num) = match typ {
                Err(e) => return Err(mlua::Error::ExternalError(Arc::new(e))),
                Ok(typ) => typ,
            };
            let n_samples = this.0.sample_count() as usize;
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
                    Err(e) => {
                        error!("info tag '{}' not found in VCF", key);
                        return Err(mlua::Error::ExternalError(Arc::new(e)));
                    }
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
                                    Ok::<LuaValue<'_>, mlua::Error>(Value::Integer(v[i]))
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
                                (_, Some(i)) => {
                                    Ok::<LuaValue<'_>, mlua::Error>(Value::Number(v[i] as f64))
                                }
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
                                        String::from_utf8_unchecked(v[i].to_vec())
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
        reg.add_method(
            "sample",
            |lua: &Lua, this: &Variant, sample_name: String| {
                let sample_id = match this.0.header().sample_id(sample_name.as_bytes()) {
                    Some(id) => id,
                    None => {
                        let msg = format!("sample '{}' not found in VCF", sample_name);
                        return Err(mlua::Error::RuntimeError(msg));
                    }
                };
                // get all format fields for this sample.
                let sample = lua.create_table().expect("error creating table");

                this.0.header().header_records().iter().for_each(|r| {
                    if let bcf::header::HeaderRecord::Format { key: _, values } = r {
                        let tag = &values["ID"];
                        let tag_bytes = tag.as_bytes();
                        let fmt = this.0.format(tag_bytes);
                        let typ = this.0.header().format_type(tag_bytes);
                        let (typ, num) = match typ {
                            Err(e) => {
                                error!("format tag '{}' error: {:?}", tag, e);
                                return;
                            }
                            Ok(typ) => typ,
                        };
                        let value = match (typ, tag_bytes) {
                            (bcf::header::TagType::Integer, _)
                            | (bcf::header::TagType::String, b"GT") => fmt
                                .integer()
                                .map(|v| match num {
                                    bcf::header::TagLength::Fixed(1) if tag_bytes != b"GT" => {
                                        Value::Integer(v[sample_id][0])
                                    }
                                    _ => {
                                        let t = lua.create_table().expect("error creating table");
                                        for (i, val) in v[sample_id].iter().enumerate() {
                                            t.raw_set(i + 1, *val).expect("error setting value");
                                        }
                                        Value::Table(t)
                                    }
                                })
                                .map_err(|e| mlua::Error::ExternalError(Arc::new(e))),
                            (bcf::header::TagType::String, _) => fmt
                                .string()
                                .map(|v| match num {
                                    bcf::header::TagLength::Fixed(1) => Value::String(
                                        lua.create_string(unsafe {
                                            String::from_utf8_unchecked(v[sample_id].to_vec())
                                        })
                                        .expect("error creating string"),
                                    ),
                                    _ => {
                                        warn!("string format tag {} is not fixed length", tag);
                                        Value::Nil
                                    }
                                })
                                .map_err(|e| mlua::Error::ExternalError(Arc::new(e))),
                            _ => Ok(Value::Nil),
                        };
                        if tag_bytes == b"GT" {
                            let gt = match value {
                                Ok(Value::Table(ref t)) => t,
                                _ => return,
                            };
                            let mut phases = vec![];
                            for i in 1..=gt.len().expect("error getting GT length") {
                                let allele = gt.get::<_, i64>(i).expect("error getting allele");
                                phases.push(allele & 1 == 1);
                                gt.raw_set(i, (allele >> 1) - 1)
                                    .expect("error setting value in GT table");
                            }
                            sample
                                .raw_set("phase", phases)
                                .expect("error setting genotype phases");
                        }
                        match value {
                            Ok(val) => sample
                                .raw_set(tag.to_string(), val)
                                .expect("error setting value"),
                            Err(e) => info!("format tag {} not found. {}", tag, e),
                        }
                    }
                });
                Ok(sample)
            },
        );
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use mlua::Lua;

    fn setup() -> (Lua, Variant) {
        let lua = Lua::new();
        register_variant(&lua).expect("error registering variant");

        let mut header = bcf::Header::new();
        header.push_record(r#"##contig=<ID=chr1,length=10000>"#.as_bytes());
        header.push_record(
            r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#.as_bytes(),
        );
        header.push_record(r#"##FILTER=<ID=PASS,Description="All filters passed">"#.as_bytes());
        header.push_record(
            r#"##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">"#.as_bytes(),
        );
        header.push_sample("NA12878".as_bytes());
        header.push_sample("NA12879".as_bytes());
        let vcf = bcf::Writer::from_path("_test.vcf", &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let _ = record.set_rid(Some(vcf.header().name2rid(b"chr1").unwrap()));
        record.set_pos(6);
        record.set_alleles(&[b"A", b"AT"]).unwrap();
        record.set_id(b"rs1234").unwrap();
        record.set_filters(&["PASS".as_bytes()]).unwrap();
        record.push_info_integer(b"DP", &[10]).unwrap();
        let alleles = &[
            bcf::record::GenotypeAllele::Unphased(0),
            bcf::record::GenotypeAllele::Phased(1),
            bcf::record::GenotypeAllele::Unphased(1),
            bcf::record::GenotypeAllele::Unphased(1),
        ];
        record.push_genotypes(alleles).unwrap();

        (lua, Variant::new(record))
    }

    #[test]
    fn test_lua_expressions() {
        let (lua, mut record) = setup();
        let globals = lua.globals();

        let expressions = vec![
            (r#"return variant.id"#, "rs1234"),
            (r#"variant.id = 'rsabc'; return variant.id"#, "rsabc"),
            (r#"return variant.REF"#, "A"),
            (r#"variant.REF = 'T'; return variant.REF"#, "T"),
            (r#"variant.ALT = {'A', 'G'}; return variant.REF"#, "T"),
            (r#"return variant.ALT[1]"#, "A"),
            (r#"return variant.ALT[2]"#, "G"),
            (r#"return variant.FILTER"#, "PASS"),
            // NOTE that we can get an integer, with 10, but we're testing
            // all strings here and verifying that the auto conversion works.
            (r#"return variant:info("DP")"#, "10"),
            // sample is 0|1 and indexing is 1-based
            (r#"s=variant:sample('NA12878'); return s.GT[1]"#, "0"),
            (r#"s=variant:sample('NA12878'); return s.GT[2]"#, "1"),
            // 2nd allele is phased to the first.
            (
                r#"s=variant:sample('NA12878'); return tostring(s.phase[2])"#,
                "true",
            ),
            // Add more expressions and expected results here...
        ];

        lua.scope(|scope| {
            let ud = scope.create_any_userdata_ref_mut(&mut record).unwrap();
            globals.raw_set("variant", ud).unwrap();

            for (expression, expected_result) in expressions {
                let exp = lua
                    .load(expression)
                    .set_name(expression)
                    .into_function()
                    .unwrap();
                let result: String = exp.call(()).unwrap();

                if result != expected_result {
                    eprintln!(
                        "expression '{}' returned '{}', expected '{}'",
                        expression, result, expected_result
                    );
                    assert_eq!(result, expected_result);
                }
            }
            Ok(())
        })
        .unwrap();
    }
}
