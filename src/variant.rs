use log::{error, info, warn};
use mlua::prelude::LuaValue;
use mlua::{AnyUserData, Lua, MetaMethod, UserDataFields, UserDataMethods, Value};
use parking_lot::Mutex;
use rust_htslib::bcf::header::{TagLength, TagType};
use rust_htslib::bcf::record::{Buffer, BufferBacked};
use rust_htslib::bcf::{self};
use rust_htslib::errors::{Error, Result};
use rustc_hash::FxHashMap;
use std::cell::RefCell;
use std::collections::HashMap;
use std::rc::Rc;
use std::sync::Arc;

/// Variant also keeps a cache of info tags to avoid repeated lookups.
pub struct HeaderMap(Rc<RefCell<FxHashMap<String, (TagType, TagLength)>>>);

impl Clone for HeaderMap {
    fn clone(&self) -> Self {
        HeaderMap(Rc::clone(&self.0))
    }
}

impl HeaderMap {
    pub fn new() -> Self {
        HeaderMap(Rc::new(RefCell::new(FxHashMap::default())))
    }
}

impl Default for HeaderMap {
    fn default() -> Self {
        HeaderMap::new()
    }
}

pub struct Variant {
    record: bcf::Record,
    header_map: HeaderMap,
}

impl Variant {
    pub fn new(record: bcf::Record, header_map: HeaderMap) -> Self {
        Variant { record, header_map }
    }
    pub fn record(&self) -> &bcf::Record {
        &self.record
    }
    pub fn header(&self) -> &bcf::header::HeaderView {
        self.record.header()
    }
    pub fn take(self) -> bcf::Record {
        self.record
    }

    pub fn info_type(&self, key: &str) -> Result<(TagType, TagLength)> {
        let t = match self.header_map.0.borrow().get(key) {
            Some((typ, num)) => return Ok((*typ, *num)),
            None => {
                let typ = self.record.header().info_type(key.as_bytes());
                match typ {
                    Err(e) => {
                        error!("info tag '{}' not found in VCF", key);
                        return Err(e);
                    }
                    Ok(t) => t,
                }
            }
        };

        self.header_map.0.borrow_mut().insert(key.to_string(), t);
        Ok(t)
    }
}

use log::{debug, log_enabled, Level};

// NEW helper functions
fn handle_format_integer<'lua>(
    lua: &'lua Lua,
    v: &BufferBacked<'_, Vec<&[i32]>, Buffer>,
    num: &bcf::header::TagLength,
    sample_id: usize,
    tag_bytes: &[u8],
) -> mlua::Result<LuaValue> {
    match num {
        bcf::header::TagLength::Fixed(1) if tag_bytes != b"GT" => {
            Ok(Value::Integer(v[sample_id][0]))
        }
        _ => {
            let t = lua
                .create_table_with_capacity(v[sample_id].len(), 0)
                .expect("error creating table");
            for (i, val) in v[sample_id].iter().enumerate() {
                t.raw_set(i + 1, *val).expect("error setting value");
            }
            Ok(Value::Table(t))
        }
    }
}

fn handle_format_float<'lua>(
    lua: &'lua Lua,
    v: &BufferBacked<'_, Vec<&[f32]>, Buffer>,
    num: &bcf::header::TagLength,
    sample_id: usize,
) -> mlua::Result<LuaValue> {
    match num {
        bcf::header::TagLength::Fixed(1) => Ok(Value::Number(v[sample_id][0] as f64)),
        _ => {
            let t = lua
                .create_table_with_capacity(v[sample_id].len(), 0)
                .expect("error creating table");
            for (i, val) in v[sample_id].iter().enumerate() {
                t.raw_set(i + 1, *val).expect("error setting value");
            }
            Ok(Value::Table(t))
        }
    }
}

fn handle_format_string<'lua>(
    lua: &'lua Lua,
    v: &BufferBacked<'_, Vec<&[u8]>, Buffer>,
    num: &bcf::header::TagLength,
    sample_id: usize,
    tag: &str,
) -> mlua::Result<LuaValue> {
    match num {
        bcf::header::TagLength::Fixed(1) => Ok(Value::String(
            lua.create_string(unsafe { String::from_utf8_unchecked(v[sample_id].to_vec()) })
                .expect("error creating string"),
        )),
        _ => {
            warn!("string format tag {} is not fixed length", tag);
            Ok(Value::Nil)
        }
    }
}

pub fn register_variant(lua: &Lua) -> mlua::Result<()> {
    lua.register_userdata_type::<Variant>(|reg| {
        reg.add_meta_function(MetaMethod::ToString, |_lua, this: AnyUserData| {
            let v = &this.borrow::<Variant>()?.record;
            let mut kstr = rust_htslib::htslib::kstring_t {
                l: 0,
                m: 0,
                s: std::ptr::null_mut(),
            };
            let h = v.header();
            unsafe { rust_htslib::htslib::vcf_format(h.inner, v.inner(), &mut kstr) };
            let s = unsafe {
                String::from_utf8_unchecked(
                    std::slice::from_raw_parts(kstr.s as *const u8, kstr.l as usize).to_vec(),
                )
            };
            eprintln!("s: {}", s);

            Ok(s)
        });
        reg.add_meta_function(
            MetaMethod::Index,
            |_lua, (_, name): (AnyUserData, String)| {
                let msg = format!("field '{}' variant.{} not found", name, name);
                Err::<LuaValue, mlua::Error>(mlua::Error::RuntimeError(msg))
            },
        );
        reg.add_field_method_get("chrom", |_, this: &Variant| {
            let c = this
                .record
                .rid()
                .map(|id| this.record.header().rid2name(id))
                .unwrap_or(Ok(b""))
                .map(|c| unsafe { String::from_utf8_unchecked(c.to_vec()) });
            c.map_err(|e| mlua::Error::ExternalError(Arc::new(e)))
        });
        reg.add_field_method_get("qual", |_, this: &Variant| Ok(this.record.qual()));
        reg.add_field_method_set("qual", |_, this: &mut Variant, val: f32| {
            this.record.set_qual(val);
            Ok(())
        });

        reg.add_field_method_get("start", |_, this: &Variant| Ok(this.record.pos()));
        reg.add_field_method_get("stop", |_, this: &Variant| Ok(this.record.end()));
        reg.add_field_method_get("pos", |_, this: &Variant| Ok(this.record.pos()));
        reg.add_field_method_set("pos", |_, this: &mut Variant, val: i64| {
            this.record.set_pos(val);
            Ok(())
        });
        reg.add_field_method_get("filters", |lua: &Lua, this: &Variant| {
            let f = this.record.filters();
            let t = lua.create_table().expect("error creating table");
            let h = this.record.header();
            for (i, id) in f.into_iter().enumerate() {
                let filter = unsafe { String::from_utf8_unchecked(h.id_to_name(id)) };
                t.raw_set(i + 1, filter).expect("error setting value");
            }
            Ok(Value::Table(t))
        });
        reg.add_field_method_set(
            "filters",
            |_, this: &mut Variant, filter: String| match this
                .record
                .set_filters(&[filter.as_bytes()])
            {
                Err(e) => Err(mlua::Error::ExternalError(Arc::new(e))),
                Ok(_) => Ok(()),
            },
        );
        reg.add_field_method_get("id", |lua: &Lua, this: &Variant| {
            let id = this.record.id();
            Ok(Value::String(unsafe {
                lua.create_string(String::from_utf8_unchecked(id.to_vec()))?
            }))
        });
        reg.add_field_method_set(
            "id",
            |_lua: &Lua, this: &mut Variant, val: String| match this.record.set_id(val.as_bytes()) {
                Err(e) => Err(mlua::Error::ExternalError(Arc::new(e))),
                Ok(_) => Ok(()),
            },
        );

        reg.add_field_method_get("REF", |lua: &Lua, this: &Variant| {
            let ref_allele = this.record.alleles()[0];
            Ok(Value::String(unsafe {
                lua.create_string(String::from_utf8_unchecked(ref_allele.to_vec()))?
            }))
        });
        reg.add_field_method_set("REF", |_lua: &Lua, this: &mut Variant, val: String| {
            let mut alleles = vec![val.as_bytes()];
            let alt_alleles = this
                .record
                .alleles()
                .iter()
                .skip(1)
                .map(|&a| a.to_owned())
                .collect::<Vec<_>>();
            alleles.extend(alt_alleles.iter().map(|a| &a[..]));

            match this.record.set_alleles(&alleles) {
                Ok(_) => Ok(()),
                Err(e) => Err(mlua::Error::ExternalError(Arc::new(e))),
            }
        });
        reg.add_field_method_set("ALT", |_lua: &Lua, this: &mut Variant, val: Vec<String>| {
            let ref_allele = this.record.alleles()[0].to_owned();
            let mut alleles = vec![&ref_allele[..]];
            alleles.extend(val.iter().map(|a| a.as_bytes()));

            match this.record.set_alleles(&alleles) {
                Ok(_) => Ok(()),
                Err(e) => Err(mlua::Error::ExternalError(Arc::new(e))),
            }
        });

        reg.add_field_method_get("ALT", |lua: &Lua, this: &Variant| {
            let alt_alleles = this.record.alleles();
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
            let f = this.record.filters();
            let h = this.record.header();
            if let Some(filter) = f.into_iter().next() {
                let filter = unsafe { String::from_utf8_unchecked(h.id_to_name(filter)) };
                return Ok(Value::String(lua.create_string(&filter)?));
            }
            Ok(Value::Nil)
        });
        reg.add_field_method_set(
            "FILTER",
            |_lua, this: &mut Variant, filter: String| match this
                .record
                .set_filters(&[filter.as_bytes()])
            {
                Err(e) => Err(mlua::Error::ExternalError(Arc::new(e))),
                Ok(_) => Ok(()),
            },
        );
        reg.add_field_method_get("genotypes", |_lua: &Lua, this: &Variant| {
            let genotypes = this.record.format(b"GT");
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
            let fmt = this.record.format(format.as_bytes());
            let typ = this.record.header().format_type(format.as_bytes());
            let (typ, num) = match typ {
                Err(e) => return Err(mlua::Error::ExternalError(Arc::new(e))),
                Ok(typ) => typ,
            };
            let n_samples = this.record.sample_count() as usize;
            let t = lua
                .create_table_with_capacity(n_samples, 0)
                .expect("error creating table");
            match typ {
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
                        Ok::<LuaValue, mlua::Error>(Value::Table(t))
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
                        Ok::<LuaValue, mlua::Error>(Value::Table(t))
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
                        Ok::<LuaValue, mlua::Error>(Value::Table(t))
                    },
                )),

                _ => unimplemented!("format type {:?}", typ),
            }
        });

        reg.add_method(
            "info",
            |lua: &Lua, this: &Variant, (key, index): (String, Option<usize>)| {
                let bkey = key.as_bytes();
                let b = Buffer::new();
                let mut info = this.record.info_shared_buffer(bkey, b);
                let typ = this.info_type(&key);
                let (typ, num) = match typ {
                    Err(e) => {
                        error!("info tag '{}' not found in VCF", key);
                        return Err(mlua::Error::ExternalError(Arc::new(e)));
                    }
                    Ok(typ) => typ,
                };
                match typ {
                    bcf::header::TagType::Integer => info
                        .integer()
                        .map(|v| handle_integer_info(lua, v, num, index))
                        .map_err(|e| mlua::Error::ExternalError(Arc::new(e))),
                    bcf::header::TagType::Float => info
                        .float()
                        .map(|v| handle_float_info(lua, v, num, index))
                        .map_err(|e| mlua::Error::ExternalError(Arc::new(e))),
                    bcf::header::TagType::String => info
                        .string()
                        .map(|v| handle_string_info(lua, v, num, index))
                        .map_err(|e| mlua::Error::ExternalError(Arc::new(e))),
                    bcf::header::TagType::Flag => info
                        .flag()
                        .map(|v| Ok::<LuaValue, mlua::Error>(Value::Boolean(v)))
                        .map_err(|e| mlua::Error::ExternalError(Arc::new(e))),
                }
            },
        );
        reg.add_method(
            "sample",
            |lua: &Lua, this: &Variant, sample_name: String| {
                let sample_id = match this.record.header().sample_id(sample_name.as_bytes()) {
                    Some(id) => id,
                    None => {
                        let msg = format!("sample '{}' not found in VCF", sample_name);
                        return Err(mlua::Error::RuntimeError(msg));
                    }
                };
                // get all format fields for this sample.
                let sample = lua.create_table().expect("error creating table");

                for r in this.record.header().header_records().iter() {
                    if let bcf::header::HeaderRecord::Format { key: _, values } = r {
                        let tag = &values["ID"];
                        let tag_bytes = tag.as_bytes();
                        let fmt = this.record.format(tag_bytes);
                        let typ = this.record.header().format_type(tag_bytes);
                        let (typ, num) = match typ {
                            Err(e) => {
                                error!("format tag '{}' error: {:?}", tag, e);
                                continue;
                            }
                            Ok(typ) => typ,
                        };

                        // Call helper functions based on TagType
                        let value = match (typ, tag_bytes) {
                            (bcf::header::TagType::Integer, _)
                            | (bcf::header::TagType::String, b"GT") => {
                                let v = fmt.integer();
                                match v {
                                    Err(Error::BcfMissingTag { tag: _, record: _ }) => {
                                        continue;
                                    }
                                    Err(e) => {
                                        error!("format tag '{}' error: {:?}", tag, e);
                                        continue;
                                    }
                                    Ok(v) => {
                                        handle_format_integer(lua, &v, &num, sample_id, tag_bytes)
                                            .map_err(|e| mlua::Error::ExternalError(Arc::new(e)))
                                    }
                                }
                            }
                            (bcf::header::TagType::Float, _) => {
                                let v = fmt.float();
                                match v {
                                    Err(Error::BcfMissingTag { tag: _, record: _ }) => {
                                        continue;
                                    }
                                    Err(e) => {
                                        error!("format tag '{}' error: {:?}", tag, e);
                                        continue;
                                    }
                                    Ok(v) => handle_format_float(lua, &v, &num, sample_id)
                                        .map_err(|e| mlua::Error::ExternalError(Arc::new(e))),
                                }
                            }
                            (bcf::header::TagType::String, _) => {
                                let v = fmt.string();
                                match v {
                                    Err(Error::BcfMissingTag { tag: _, record: _ }) => continue,
                                    Err(e) => {
                                        error!("format tag '{}' error: {:?}", tag, e);
                                        continue;
                                    }
                                    Ok(v) => handle_format_string(lua, &v, &num, sample_id, tag)
                                        .map_err(|e| mlua::Error::ExternalError(Arc::new(e))),
                                }
                            }

                            _ => Ok(Value::Nil),
                        };

                        if tag_bytes == b"GT" {
                            let gt = match value {
                                Ok(Value::Table(ref t)) => t,
                                _ => continue,
                            };
                            let mut phases = vec![];
                            let mut alts = 0;
                            for i in 1..=gt.len().expect("error getting GT length") {
                                let allele = gt.get::<i64>(i).expect("error getting allele");
                                phases.push(allele & 1 == 1);
                                alts += (allele >> 1) - 1;
                                gt.raw_set(i, (allele >> 1) - 1)
                                    .expect("error setting value in GT table");
                            }
                            sample
                                .raw_set("phase", phases)
                                .expect("error setting genotype phases");
                            sample
                                .raw_set("alts", alts)
                                .expect("error setting genotype alts");
                        }
                        match value {
                            Ok(val) => sample
                                .raw_set(tag.to_string(), val)
                                .expect("error setting value"),
                            Err(e) => info!("format tag {} not found. {}", tag, e),
                        }
                    }
                }
                Ok(sample)
            },
        );

        reg.add_method(
            "samples",
            |lua: &Lua, this: &Variant, fields: Option<HashMap<String, bool>>| {
                let mut samples = HashMap::new();

                // Create a table for each sample
                let sample_names = this
                    .record
                    .header()
                    .samples()
                    .iter()
                    .map(|s| unsafe { String::from_utf8_unchecked(s.to_vec()) })
                    .collect::<Vec<_>>();
                for sample_name in sample_names.iter() {
                    samples.insert(sample_name.to_string(), HashMap::new());
                }

                // Process all format fields
                for r in this.record.header().header_records().iter() {
                    if let bcf::header::HeaderRecord::Format { key: _, values } = r {
                        let tag = &values["ID"];
                        let tag_bytes = tag.as_bytes();

                        // Skip if not GT and not in requested fields
                        if tag_bytes != b"GT" {
                            if let Some(ref fields) = fields {
                                let should_include =
                                    fields.get(tag.to_string().as_str()).unwrap_or(&false);
                                if !should_include {
                                    continue;
                                }
                            }
                        }

                        let fmt = this.record.format(tag_bytes);
                        let typ = this.record.header().format_type(tag_bytes);
                        let (typ, num) = match typ {
                            Err(e) => {
                                error!("format tag '{}' error: {:?}", tag, e);
                                continue;
                            }
                            Ok(typ) => typ,
                        };

                        // Process each format field type
                        match (typ, tag_bytes) {
                            (bcf::header::TagType::Integer, _)
                            | (bcf::header::TagType::String, b"GT") => {
                                let v = fmt.integer();
                                match v {
                                    Err(Error::BcfMissingTag { tag: _, record: _ }) => continue,
                                    Err(e) => {
                                        error!("format tag '{}' error: {:?}", tag, e);
                                        continue;
                                    }
                                    Ok(v) => {
                                        for (sample_id, sample_name) in
                                            sample_names.iter().enumerate()
                                        {
                                            let sample = samples
                                                .get_mut(sample_name)
                                                .expect("error getting sample map");

                                            let value = handle_format_integer(
                                                lua, &v, &num, sample_id, tag_bytes,
                                            )
                                            .expect("error handling integer format");

                                            if tag_bytes == b"GT" {
                                                if let Value::Table(ref gt) = value {
                                                    let mut phases = Vec::with_capacity(2);
                                                    let mut alts = 0;
                                                    for i in 1..=gt
                                                        .len()
                                                        .expect("error getting GT length")
                                                    {
                                                        let allele = gt
                                                            .get::<i64>(i)
                                                            .expect("error getting allele");
                                                        phases.push(allele & 1 == 1);
                                                        alts += (allele >> 1) - 1;
                                                        gt.raw_set(i, (allele >> 1) - 1).expect(
                                                            "error setting value in GT table",
                                                        );
                                                    }
                                                    let phases_table = lua
                                                        .create_table()
                                                        .expect("error creating table");
                                                    for (i, phase) in phases.iter().enumerate() {
                                                        phases_table
                                                            .raw_set(i + 1, *phase)
                                                            .expect("error setting phase");
                                                    }
                                                    sample.insert(
                                                        "phase".to_string(),
                                                        Value::Table(phases_table),
                                                    );
                                                    sample.insert(
                                                        "alts".to_string(),
                                                        Value::Integer(alts as i32),
                                                    );
                                                }
                                            }
                                            sample.insert(tag.to_string(), value);
                                        }
                                    }
                                }
                            }
                            (bcf::header::TagType::Float, _) => {
                                let v = fmt.float();
                                match v {
                                    Err(Error::BcfMissingTag { tag: _, record: _ }) => continue,
                                    Err(e) => {
                                        error!("format tag '{}' error: {:?}", tag, e);
                                        continue;
                                    }
                                    Ok(v) => {
                                        for (sample_id, sample_name) in
                                            sample_names.iter().enumerate()
                                        {
                                            let sample = samples
                                                .get_mut(sample_name)
                                                .expect("error getting sample map");

                                            let value =
                                                handle_format_float(lua, &v, &num, sample_id)
                                                    .expect("error handling float format");
                                            sample.insert(tag.to_string(), value);
                                        }
                                    }
                                }
                            }
                            (bcf::header::TagType::String, _) => {
                                let v = fmt.string();
                                match v {
                                    Err(Error::BcfMissingTag { tag: _, record: _ }) => continue,
                                    Err(e) => {
                                        error!("format tag '{}' error: {:?}", tag, e);
                                        continue;
                                    }
                                    Ok(v) => {
                                        for (sample_id, sample_name) in
                                            sample_names.iter().enumerate()
                                        {
                                            let sample = samples
                                                .get_mut(sample_name)
                                                .expect("error getting sample map");

                                            let value =
                                                handle_format_string(lua, &v, &num, sample_id, tag)
                                                    .expect("error handling string format");
                                            sample.insert(tag.to_string(), value);
                                        }
                                    }
                                }
                            }
                            _ => continue,
                        }
                    }
                }
                Ok(samples)
            },
        );
    })
}

fn handle_integer_info<'lua>(
    lua: &'lua Lua,
    v: Option<BufferBacked<'_, &[i32], Buffer>>,
    num: TagLength,
    index: Option<usize>,
) -> mlua::Result<LuaValue> {
    match v {
        Some(v) => match (num, index) {
            (bcf::header::TagLength::Fixed(1), None) => Ok(Value::Integer(v[0])),
            (_, Some(i)) => Ok(Value::Integer(v[i])),
            _ => {
                let t = lua.create_table()?;
                for (i, val) in v.iter().enumerate() {
                    t.raw_set(i + 1, *val)?;
                }
                Ok(Value::Table(t))
            }
        },
        None => Ok(Value::Nil),
    }
}

fn handle_float_info<'lua>(
    lua: &'lua Lua,
    v: Option<BufferBacked<'_, &[f32], Buffer>>,
    num: TagLength,
    index: Option<usize>,
) -> mlua::Result<LuaValue> {
    match v {
        Some(v) => match (num, index) {
            (bcf::header::TagLength::Fixed(1), None) => Ok(Value::Number(f64::from(v[0]))),
            (_, Some(i)) => Ok(Value::Number(f64::from(v[i]))),
            _ => {
                let t = lua.create_table()?;
                for (i, val) in v.iter().enumerate() {
                    t.raw_set(i + 1, f64::from(*val))?;
                }
                Ok(Value::Table(t))
            }
        },
        None => Ok(Value::Nil),
    }
}

fn handle_string_info<'lua>(
    lua: &'lua Lua,
    v: Option<BufferBacked<'_, Vec<&[u8]>, Buffer>>,
    num: TagLength,
    index: Option<usize>,
) -> mlua::Result<LuaValue> {
    match v {
        Some(v) => match (num, index) {
            (bcf::header::TagLength::Fixed(1), None) => {
                Ok(Value::String(lua.create_string(unsafe {
                    String::from_utf8_unchecked(v[0].to_vec())
                })?))
            }
            (_, Some(i)) => {
                Ok(Value::String(lua.create_string(unsafe {
                    String::from_utf8_unchecked(v[i].to_vec())
                })?))
            }
            _ => {
                let t = lua.create_table()?;
                for (i, s) in v.iter().enumerate() {
                    t.raw_set(i + 1, unsafe { String::from_utf8_unchecked(s.to_vec()) })?;
                }
                Ok(Value::Table(t))
            }
        },
        None => Ok(Value::Nil),
    }
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
        // Add Format fields for testing
        header.push_record(
            r#"##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">"#.as_bytes(),
        );
        header.push_record(
            r#"##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality">"#.as_bytes(),
        );
        header.push_record(
            r#"##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">"#.as_bytes(),
        );
        header.push_record(
            r#"##FORMAT=<ID=SQ,Number=1,Type=String,Description="String Quality">"#.as_bytes(),
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

        // Push sample-specific format data.
        record.push_format_integer(b"DP", &[11, 12]).unwrap();
        record.push_format_float(b"GQ", &[40.0, 50.0]).unwrap();
        record
            .push_format_integer(b"HQ", &[10, 20, 30, 40])
            .unwrap(); // 2 values per sample
                       //record.push_format_string(b"SQ", &[b"abc", b"def"]).unwrap();

        (lua, Variant::new(record, HeaderMap::new()))
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
            // Test Integer Format Field
            (r#"s=variant:sample('NA12878'); return s.DP"#, "11"),
            (r#"s=variant:sample('NA12879'); return s.DP"#, "12"),
            // Test Float Format Field
            (
                r#"s=variant:sample('NA12878'); return s.GQ"#,
                "40", // Lua converts to string
            ),
            (r#"s=variant:sample('NA12879'); return s.GQ"#, "50"),
            // Test String Format Field
            //(r#"s=variant:sample('NA12878'); return s.SQ"#, "abc"),
            //(r#"s=variant:sample('NA12879'); return s.SQ"#, "def"),
            // Test multi-value Integer Format Field
            (r#"s=variant:sample('NA12878'); return s.HQ[1]"#, "10"),
            (r#"s=variant:sample('NA12878'); return s.HQ[2]"#, "20"),
            (r#"s=variant:sample('NA12879'); return s.HQ[1]"#, "30"),
            (r#"s=variant:sample('NA12879'); return s.HQ[2]"#, "40"),
            // Test samples() method
            (r#"s=variant:samples(); return s.NA12878.GT[1]"#, "0"),
            (r#"s=variant:samples(); return s.NA12878.GT[2]"#, "1"),
            (r#"s=variant:samples(); return s.NA12879.GT[1]"#, "1"),
            (r#"s=variant:samples(); return s.NA12879.GT[2]"#, "1"),
            (
                r#"s=variant:samples(); return tostring(s.NA12878.phase[2])"#,
                "true",
            ),
            (r#"s=variant:samples(); return s.NA12878.DP"#, "11"),
            (r#"s=variant:samples(); return s.NA12879.DP"#, "12"),
            (r#"s=variant:samples(); return s.NA12878.GQ"#, "40"),
            (r#"s=variant:samples(); return s.NA12879.GQ"#, "50"),
            (r#"s=variant:samples(); return s.NA12878.HQ[1]"#, "10"),
            (r#"s=variant:samples(); return s.NA12878.HQ[2]"#, "20"),
            (r#"s=variant:samples(); return s.NA12879.HQ[1]"#, "30"),
            (r#"s=variant:samples(); return s.NA12879.HQ[2]"#, "40"),
            // Test samples() method with field filtering
            (r#"s=variant:samples({DP=true}); return s.NA12878.DP"#, "11"),
            (r#"s=variant:samples({DP=true}); return s.NA12879.DP"#, "12"),
            // Test that GT is always included even when not specified
            (
                r#"s=variant:samples({DP=true}); return s.NA12878.GT[1]"#,
                "0",
            ),
            (
                r#"s=variant:samples({DP=true}); return s.NA12878.GT[2]"#,
                "1",
            ),
            // Test that non-requested fields are not included
            (
                r#"s=variant:samples({DP=true}); return tostring(s.NA12878.GQ)"#,
                "nil",
            ),
            (
                r#"s=variant:samples({DP=true}); return tostring(s.NA12878.HQ)"#,
                "nil",
            ),
            // Test multiple requested fields
            (
                r#"s=variant:samples({DP=true,GQ=true}); return s.NA12878.GQ"#,
                "40",
            ),
            (
                r#"s=variant:samples({DP=true,GQ=true}); return s.NA12879.GQ"#,
                "50",
            ),
            // Test that phase and alts are included with GT
            (
                r#"s=variant:samples({DP=true}); return tostring(s.NA12878.phase[2])"#,
                "true",
            ),
            (
                r#"s=variant:samples({DP=true}); return s.NA12878.alts"#,
                "1",
            ),
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
                let result: String = exp
                    .call(())
                    .expect(&format!("error calling expression: {}", expression));

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
