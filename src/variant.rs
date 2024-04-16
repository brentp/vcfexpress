use mlua::prelude::LuaValue;
use mlua::{AnyUserData, Lua, MetaMethod, UserDataFields, UserDataMethods, Value};
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
