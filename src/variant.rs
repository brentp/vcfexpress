use log::{error, info, warn};
use rusty_v8 as v8;
use rust_htslib::bcf::header::{TagLength, TagType};
use rust_htslib::bcf::record::Buffer;
use rust_htslib::bcf::{self};
use rust_htslib::errors::Result;
use rustc_hash::FxHashMap;
use std::cell::RefCell;
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

// Helper function to wrap Variant as an internal field
fn create_variant_object<'a>(
    scope: &mut v8::HandleScope<'a>,
    variant: Arc<Variant>,
) -> v8::Local<'a, v8::Object> {
    // Create an object template with one internal field
    let object_template = v8::ObjectTemplate::new(scope);
    object_template.set_internal_field_count(1);

    // Create an instance of the template
    let object = object_template.new_instance(scope).unwrap();

    // Create an external reference to your Rust struct
    let external_variant = v8::External::new(scope, Arc::into_raw(variant) as *mut _);

    // Store the external reference in the internal field
    object.set_internal_field(0, external_variant.into());

    object
}

fn get_start(
    scope: &mut v8::HandleScope,
    _name: v8::Local<v8::Name>,
    args: v8::PropertyCallbackArguments,
    mut rv: v8::ReturnValue,
) {
    let this = args.this();

    // Get the Variant from the internal field
    let internal_field = this.get_internal_field(scope, 0).unwrap();
    let external_variant = v8::Local::<v8::External>::try_from(internal_field).unwrap();
    let variant = unsafe { &*(external_variant.value() as *const Variant) };

    // Return the `start` field value to JavaScript
    rv.set(v8::Integer::new(scope, variant.record.pos() as i32).into());
}


pub (crate) fn register_variant<'a>(scope: &mut v8::HandleScope<'a>, object: v8::Local<'a, v8::Object>) {
    let start_name = v8::String::new(scope, "start").unwrap();

    // Define the property with getter and setter for `start`
    object.set_accessor(
        scope,
        start_name.into(),
        get_start,
        None, //Some(start_setter),
    );
}



/*
// Implement methods
fn info_method(
    scope: &mut v8::HandleScope,
    args: v8::FunctionCallbackArguments,
    mut retval: v8::ReturnValue,
) {
    let variant = args.this().get_internal_field(scope, 0).unwrap().is_external().unwrap();
    let variant = unsafe { &*(variant.value() as *const Variant) };
    
    if args.length() < 1 {
        return;
    }
    
    let key = args.get(0).to_string(scope).unwrap().to_rust_string_lossy(scope);
    
    // Implement info retrieval logic here
    // ...

    // Set the return value based on the info type
    // retval.set(...);
}
    */

// ... implement other methods

#[cfg(test)]
mod tests {
    use super::*;
    use v8;

    fn setup() -> (v8::OwnedIsolate, Variant) {
        let mut isolate = v8::Isolate::new(v8::CreateParams::default());
        let context = {
            let handle_scope = &mut v8::HandleScope::new(&mut isolate);
            v8::Context::new(handle_scope)
        };
        let scope = &mut v8::HandleScope::with_context(&mut isolate, &context);

        register_variant(&mut isolate, &context).expect("error registering variant");

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
        
        let variant = Variant::new(record, HeaderMap::new());

        (isolate, variant)
    }

    #[test]
    fn test_javascript_expressions() {
        let (mut isolate, record) = setup();
        let context = {
            let handle_scope = &mut v8::HandleScope::new(&mut isolate);
            v8::Context::new(handle_scope)
        };
        let scope = &mut v8::HandleScope::with_context(&mut isolate, &context);
        let global = context.global(scope);

        let expressions = vec![
            (r#"return variant.start"#, "6"),
            /*
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
            */
        ];

        let ud = v8::External::new(scope, &record as *const Variant as *mut std::ffi::c_void);
        global.set(
            scope,
            v8::String::new(scope, "variant").unwrap().into(),
            ud.into(),
        ).unwrap();

        for (expression, expected_result) in expressions {
            let script = v8::Script::compile(scope, v8::String::new(scope, expression).unwrap(), None).unwrap();
            let result = script.run(scope).unwrap();
            let result = result.to_string(scope).unwrap().to_rust_string_lossy(scope);

            if result != expected_result {
                eprintln!(
                    "expression '{}' returned '{}', expected '{}'",
                    expression, result, expected_result
                );
                assert_eq!(result, expected_result);
            }
        }
    }
}
