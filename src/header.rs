use rusty_v8 as v8;
use rust_htslib::bcf::header::{Header, HeaderView};
use rust_htslib::bcf::HeaderRecord;
use std::collections::HashMap;

fn handle_hash_get<'a>(
    tbl: &'a HashMap<String, String>,
    key: &str,
    func: &str,
) -> Result<&'a str, Box<dyn std::error::Error>> {
    match tbl.get(key) {
        Some(x) => Ok(x),
        None => Err(Box::new(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            format!(
                "must specify {} in the argument to add_{}. got {:?}",
                key, func, tbl
            ),
        ))),
    }
}

macro_rules! find_record_match {
    ($record:expr, $key:expr, $record_type:ident) => {
        match $record {
            HeaderRecord::$record_type { key: _, values } => {
                if values.get("ID") != Some(&$key.to_string()) {
                    return None;
                }
                Some(HashMap::from_iter(
                    values
                        .into_iter()
                        .map(|(k, v)| (k.to_string(), v.to_string())),
                ))
            }
            _ => None,
        }
    };
}

fn find_record(
    records: &[HeaderRecord],
    key: &str,
    hdr_type: ::libc::c_uint,
) -> Result<HashMap<String, String>, Box<dyn std::error::Error>> {
    let hrec = records
        .iter()
        .filter_map(|x| {
            if hdr_type == rust_htslib::htslib::BCF_HL_INFO {
                find_record_match!(x, key, Info)
            } else if hdr_type == rust_htslib::htslib::BCF_HL_FMT {
                find_record_match!(x, key, Format)
            } else {
                None
            }
        })
        .next();

    if let Some(hrec) = hrec {
        return Ok(hrec);
    }
    Err(Box::new(std::io::Error::new(
        std::io::ErrorKind::NotFound,
        format!("key {}, hdr_type:{:?} not found in header", key, hdr_type),
    )))
}

pub(crate) fn register_header(isolate: &mut v8::Isolate, context: &v8::Local<v8::Context>) -> Result<(), Box<dyn std::error::Error>> {
    let scope = &mut v8::HandleScope::with_context(isolate, context);
    let global = context.global(scope);

    let header_template = v8::ObjectTemplate::new(scope);
    
    // Add methods
    header_template.set(
        v8::String::new(scope, "info_get").unwrap().into(),
        v8::FunctionTemplate::new(scope, info_get_method).into(),
    );
    header_template.set(
        v8::String::new(scope, "format_get").unwrap().into(),
        v8::FunctionTemplate::new(scope, format_get_method).into(),
    );
    header_template.set(
        v8::String::new(scope, "add_info").unwrap().into(),
        v8::FunctionTemplate::new(scope, add_info_method).into(),
    );
    header_template.set(
        v8::String::new(scope, "add_format").unwrap().into(),
        v8::FunctionTemplate::new(scope, add_format_method).into(),
    );
    // ... add more methods

    global.set(
        scope,
        v8::String::new(scope, "Header").unwrap().into(),
        header_template.into(),
    ).unwrap();

    Ok(())
}

fn info_get_method(
    scope: &mut v8::HandleScope,
    args: v8::FunctionCallbackArguments,
    mut retval: v8::ReturnValue,
) {
    let header = args.this().get_internal_field(scope, 0).unwrap().as_external().unwrap();
    let header = unsafe { &*(header.value() as *const HeaderView) };
    
    if args.length() < 1 {
        return;
    }
    
    let find_key = args.get(0).to_string(scope).unwrap().to_rust_string_lossy(scope);
    
    // Implement info_get logic here
    // ...

    // Set the return value
    // retval.set(...);
}

fn format_get_method(
    scope: &mut v8::HandleScope,
    args: v8::FunctionCallbackArguments,
    mut retval: v8::ReturnValue,
) {
    let header = args.this().get_internal_field(scope, 0).unwrap().as_external().unwrap();
    let header = unsafe { &*(header.value() as *const HeaderView) };
    
    if args.length() < 1 {
        return;
    }
    
    let find_key = args.get(0).to_string(scope).unwrap().to_rust_string_lossy(scope);
    
    // Implement format_get logic here
    // ...

    // Set the return value
    // retval.set(...);
}

fn add_info_method(
    scope: &mut v8::HandleScope,
    args: v8::FunctionCallbackArguments,
    mut retval: v8::ReturnValue,
) {
    let header = args.this().get_internal_field(scope, 0).unwrap().as_external().unwrap();
    let header = unsafe { &mut *(header.value() as *mut HeaderView) };
    
    if args.length() < 1 {
        return;
    }
    
    let tbl = args.get(0).to_object(scope).unwrap();
    let mut tbl_map = HashMap::new();
    let keys = tbl.get_own_property_names(scope).unwrap();
    for i in 0..keys.length() {
        let key = keys.get(scope, i).unwrap().to_string(scope).unwrap().to_rust_string_lossy(scope);
        let value = tbl.get(scope, v8::String::new(scope, &key).unwrap().into()).unwrap().to_string(scope).unwrap().to_rust_string_lossy(scope);
        tbl_map.insert(key, value);
    }
    
    let c_str = std::ffi::CString::new(format!(
        r#"##INFO=<ID={},Number={},Type={},Description={}>"#,
        handle_hash_get(&tbl_map, "ID", "info").unwrap(),
        handle_hash_get(&tbl_map, "Number", "info").unwrap(),
        handle_hash_get(&tbl_map, "Type", "info").unwrap(),
        handle_hash_get(&tbl_map, "Description", "info").unwrap(),
    ))
    .expect("CString::new failed");
    let ret =
        unsafe { rust_htslib::htslib::bcf_hdr_append(header.inner, c_str.as_ptr()) };
    if ret != 0 {
        log::error!("Error adding INFO field for {:?}: {}", tbl_map, ret);
        return;
    }
    let ret = unsafe { rust_htslib::htslib::bcf_hdr_sync(header.inner) };
    if ret != 0 {
        log::warn!(
            "Error syncing header after adding INFO field for {:?}: {}",
            tbl_map,
            ret
        );
    }
}

fn add_format_method(
    scope: &mut v8::HandleScope,
    args: v8::FunctionCallbackArguments,
    mut retval: v8::ReturnValue,
) {
    let header = args.this().get_internal_field(scope, 0).unwrap().as_external().unwrap();
    let header = unsafe { &mut *(header.value() as *mut HeaderView) };
    
    if args.length() < 1 {
        return;
    }
    
    let tbl = args.get(0).to_object(scope).unwrap();
    let mut tbl_map = HashMap::new();
    let keys = tbl.get_own_property_names(scope).unwrap();
    for i in 0..keys.length() {
        let key = keys.get(scope, i).unwrap().to_string(scope).unwrap().to_rust_string_lossy(scope);
        let value = tbl.get(scope, v8::String::new(scope, &key).unwrap().into()).unwrap().to_string(scope).unwrap().to_rust_string_lossy(scope);
        tbl_map.insert(key, value);
    }
    
    let c_str = std::ffi::CString::new(format!(
        r#"##FORMAT=<ID={},Number={},Type={},Description="{}">"#,
        handle_hash_get(&tbl_map, "ID", "format").unwrap(),
        handle_hash_get(&tbl_map, "Number", "format").unwrap(),
        handle_hash_get(&tbl_map, "Type", "format").unwrap(),
        handle_hash_get(&tbl_map, "Description", "format").unwrap(),
    ))
    .expect("CString::new failed");
    let ret =
        unsafe { rust_htslib::htslib::bcf_hdr_append(header.inner, c_str.as_ptr()) };
    if ret != 0 {
        log::error!("Error adding FORMAT field for {:?}: {}", tbl_map, ret);
        return;
    }
    _ = unsafe { rust_htslib::htslib::bcf_hdr_sync(header.inner) };
}

#[cfg(test)]
mod tests {
    use super::*;
    use rusty_v8 as v8;

    fn setup() -> (v8::OwnedIsolate, Header, HeaderView) {
        let mut isolate = v8::Isolate::new(v8::CreateParams::default());
        let handle_scope = &mut v8::HandleScope::new(&mut isolate);
        let context = v8::Context::new(handle_scope);

        register_header(&mut isolate, &context).unwrap();

        let mut header = Header::new();
        header.push_record(r#"##contig=<ID=chr1,length=10000>"#.as_bytes());

        header.push_sample("Sample1".as_bytes());
        header.push_sample("Sample2".as_bytes());
        unsafe { rust_htslib::htslib::bcf_hdr_sync(header.inner) };
        let header_t = unsafe { rust_htslib::htslib::bcf_hdr_dup(header.inner) };
        let header_view = HeaderView::new(header_t);

        (isolate, header, header_view)
    }

    #[test]
    fn test_v8_header_samples() {
        let (isolate, _header, mut header_view) = setup();
        let scope = &mut v8::HandleScope::new(&isolate);
        let context = v8::Context::new(scope);

        let exp = v8::Script::compile(
            scope,
            r#"
            // TODO: this is broken
            header.samples = ["Sample1"];
            return header.samples.join(",");
            "#,
            None,
        )
        .unwrap();

        let result = exp.run(scope).unwrap().unwrap().to_string(scope).unwrap().to_rust_string_lossy(scope);
        assert_eq!(result, "Sample1");
    }

    #[test]
    fn test_add_info() {
        let (isolate, _header, mut header_view) = setup();
        let scope = &mut v8::HandleScope::new(&isolate);
        let context = v8::Context::new(scope);

        let exp = v8::Script::compile(
            scope,
            r#"
            header.add_info({ID:"TEST", Number:"1", Type:"Integer", Description:"Test field"});
            return header.toString();
            "#,
            None,
        )
        .unwrap();

        let result = exp.run(scope).unwrap().unwrap().to_string(scope).unwrap().to_rust_string_lossy(scope);
        let expected = "##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description=\"All filters passed\">\n##contig=<ID=chr1,length=10000>\n##INFO=<ID=TEST,Number=1,Type=Integer,Description=Test field>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2\n";
        assert_eq!(
            result,
            expected,
        );
    }

    #[test]
    fn test_add_format() {
        let (isolate, _header, mut header_view) = setup();
        let scope = &mut v8::HandleScope::new(&isolate);
        let context = v8::Context::new(scope);

        let exp = v8::Script::compile(
            scope,
            r#"
            header.add_format({ID:"TEST", Number:"1", Type:"Integer", Description:"Test field"});
            return header.toString();
            "#,
            None,
        )
        .unwrap();

        let result = exp.run(scope).unwrap().unwrap().to_string(scope).unwrap().to_rust_string_lossy(scope);
        let expected = "##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description=\"All filters passed\">\n##contig=<ID=chr1,length=10000>\n##FORMAT=<ID=TEST,Number=1,Type=Integer,Description=\"Test field\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2\n";
        assert_eq!(
            result,
            expected,
        );
    }
}
