use mlua::{AnyUserData, Lua, MetaMethod, UserDataFields, UserDataMethods};
use rust_htslib::bcf::header::{Header, HeaderView};
use std::collections::HashMap;
use std::sync::Arc;

fn handle_hash_get<'a>(
    tbl: &'a HashMap<String, String>,
    key: &str,
    func: &str,
) -> Result<&'a str, mlua::Error> {
    match tbl.get(key) {
        Some(x) => Ok(x),
        None => {
            Err(mlua::Error::ExternalError(Arc::new(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!(
                    "must specify {} in the argument to add_{}. got {:?}",
                    key, func, tbl
                ),
            ))))
        }
    }
}

pub(crate) fn register_header(lua: &Lua) -> mlua::Result<()> {
    lua.register_userdata_type::<HeaderView>(|reg| {
        reg.add_meta_function(MetaMethod::ToString, |_lua, this: AnyUserData| {
            let this = this.borrow::<HeaderView>()?;

            let mut kstr = rust_htslib::htslib::kstring_t {
                l: 0,
                m: 0,
                s: std::ptr::null_mut(),
            };
            if unsafe { rust_htslib::htslib::bcf_hdr_format(this.inner, 0, &mut kstr) } != 0 {
                return Err(mlua::Error::ExternalError(Arc::new(
                    std::io::Error::last_os_error(),
                )));
            }
            let s = unsafe {
                String::from_utf8_unchecked(
                    std::slice::from_raw_parts(kstr.s as *const u8, kstr.l as usize).to_vec(),
                )
            };

            Ok(s)
        });
        reg.add_field_method_get("samples", |_lua, this: &HeaderView| {
            let samples = this
                .samples()
                .iter()
                .map(|&x| String::from_utf8_lossy(x).to_string())
                .collect::<Vec<_>>();
            Ok(samples)
        });
        reg.add_field_method_set(
            "samples",
            |_lua, this: &mut HeaderView, samples: Vec<String>| {
                let sample_bytes = samples.iter().map(|x| x.as_bytes()).collect::<Vec<_>>();
                match Header::from_template_subset(this, &sample_bytes) {
                    Ok(h) => {
                        //_ = unsafe { rust_htslib::htslib::bcf_hdr_sync(h.inner) };
                        let header_t = unsafe { rust_htslib::htslib::bcf_hdr_dup(h.inner) };
                        *this = HeaderView::new(header_t);
                        eprintln!(
                            "samples from c directly: {:?}",
                            this.samples()
                                .iter()
                                .map(|&x| String::from_utf8_lossy(x).to_string())
                                .collect::<Vec<_>>()
                        );
                        Ok(())
                    }
                    Err(e) => Err(mlua::Error::ExternalError(Arc::new(e))),
                }
            },
        );
        reg.add_function_mut(
            "add_info",
            |_lua, (ud, tbl): (AnyUserData, HashMap<String, String>)| {
                let this = ud.borrow_mut::<HeaderView>()?;
                let c_str = std::ffi::CString::new(format!(
                    r#"##INFO=<ID={},Number={},Type={},Description="{}">"#,
                    handle_hash_get(&tbl, "ID", "info")?,
                    handle_hash_get(&tbl, "Number", "info")?,
                    handle_hash_get(&tbl, "Type", "info")?,
                    handle_hash_get(&tbl, "Description", "info")?,
                ))
                .expect("CString::new failed");
                let ret =
                    unsafe { rust_htslib::htslib::bcf_hdr_append(this.inner, c_str.as_ptr()) };
                if ret != 0 {
                    log::error!("Error adding INFO field for {:?}: {}", tbl, ret);
                    return Err(mlua::Error::ExternalError(Arc::new(
                        std::io::Error::last_os_error(),
                    )));
                }
                _ = unsafe { rust_htslib::htslib::bcf_hdr_sync(this.inner) };
                Ok(())
            },
        );
        reg.add_function_mut(
            "add_format",
            |_lua, (ud, tbl): (AnyUserData, HashMap<String, String>)| {
                let this = ud.borrow_mut::<HeaderView>()?;
                let c_str = std::ffi::CString::new(format!(
                    r#"##FORMAT=<ID={},Number={},Type={},Description="{}">"#,
                    handle_hash_get(&tbl, "ID", "format")?,
                    handle_hash_get(&tbl, "Number", "format")?,
                    handle_hash_get(&tbl, "Type", "format")?,
                    handle_hash_get(&tbl, "Description", "format")?,
                ))
                .expect("CString::new failed");
                let ret =
                    unsafe { rust_htslib::htslib::bcf_hdr_append(this.inner, c_str.as_ptr()) };
                if ret != 0 {
                    log::error!("Error adding FORMAT field for {:?}: {}", tbl, ret);
                    return Err(mlua::Error::ExternalError(Arc::new(
                        std::io::Error::last_os_error(),
                    )));
                }
                _ = unsafe { rust_htslib::htslib::bcf_hdr_sync(this.inner) };
                Ok(())
            },
        );
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use mlua::Lua;

    fn setup() -> (Lua, Header, HeaderView) {
        let lua = Lua::new();
        register_header(&lua).unwrap();

        let mut header = Header::new();
        header.push_record(r#"##contig=<ID=chr1,length=10000>"#.as_bytes());

        header.push_sample("Sample1".as_bytes());
        header.push_sample("Sample2".as_bytes());
        unsafe { rust_htslib::htslib::bcf_hdr_sync(header.inner) };
        let header_t = unsafe { rust_htslib::htslib::bcf_hdr_dup(header.inner) };
        let header_view = HeaderView::new(header_t);

        (lua, header, header_view)
    }

    #[test]
    fn test_lua_header_samples() {
        let (lua, _header, mut header_view) = setup();
        let globals = lua.globals();

        let exp = lua
            .load(
                r#"
            -- TODO: this is broken
            header.samples = {"Sample1"};
            return table.concat(header.samples, ",")
            "#,
            )
            .set_name("test_lua_header_samples")
            .into_function()
            .expect("error in test_lua_header_samples");

        lua.scope(|scope| {
            globals.set(
                "header",
                scope.create_any_userdata_ref_mut(&mut header_view)?,
            )?;
            let result: String = exp.call(())?;
            assert_eq!(result, "Sample1");

            Ok(())
        })
        .expect("error in test_lua_header_samples")
    }

    #[test]
    fn test_add_info() {
        let (lua, _header, mut header_view) = setup();
        let globals = lua.globals();

        let exp = lua
            .load(
                r#"
            header:add_info({ID="TEST", Number="1", Type="Integer", Description="Test field"});
            return tostring(header)
            "#,
            )
            .set_name("test_add_info")
            .into_function()
            .expect("error in test_add_info");

        lua.scope(|scope| {
            globals.set(
                "header",
                scope.create_any_userdata_ref_mut(&mut header_view)?,
            )?;
            let result: String = exp.call(())?;
            let expected = "##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description=\"All filters passed\">\n##contig=<ID=chr1,length=10000>\n##INFO=<ID=TEST,Number=1,Type=Integer,Description=\"Test field\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2\n";
            assert_eq!(
                result,
                expected,
            );

            Ok(())
        })
        .expect("error in test_add_info")
    }

    #[test]
    fn test_add_format() {
        let (lua, _header, mut header_view) = setup();
        let globals = lua.globals();

        let exp = lua
            .load(
                r#"
            header:add_format({ID="TEST", Number="1", Type="Integer", Description="Test field"});
            return tostring(header)
            "#,
            )
            .set_name("test_add_format")
            .into_function()
            .expect("error in test_add_format");

        lua.scope(|scope| {
            globals.set(
                "header",
                scope.create_any_userdata_ref_mut(&mut header_view)?,
            )?;
            let result: String = exp.call(())?;
            let expected = "##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description=\"All filters passed\">\n##contig=<ID=chr1,length=10000>\n##FORMAT=<ID=TEST,Number=1,Type=Integer,Description=\"Test field\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2\n";
            assert_eq!(
                result,
                expected,
            );

            Ok(())
        })
        .expect("error in test_add_format")
    }
}
