use mlua::{AnyUserData, Lua, MetaMethod, UserDataFields, UserDataMethods};
use rust_htslib::bcf::header::{Header, HeaderView};
use std::sync::Arc;

pub(crate) fn register_header(lua: &Lua) -> mlua::Result<()> {
    lua.register_userdata_type::<HeaderView>(|reg| {
        reg.add_meta_function(MetaMethod::ToString, |_lua, this: AnyUserData| {
            let this = this.borrow::<HeaderView>()?;
            Ok(format!("Header: {:?}", this))
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
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use mlua::Lua;

    #[test]
    fn test_lua_header_samples() {
        let lua = Lua::new();
        let globals = lua.globals();
        register_header(&lua).unwrap();
        //let header_view = create_sample_header();
        let mut header = Header::new();
        header.push_record(r#"##contig=<ID=chr1,length=10000>"#.as_bytes());

        header.push_sample("Sample1".as_bytes());
        header.push_sample("Sample2".as_bytes());
        // this is required because rust-htslib doesn't call bcf_hdr_sync
        //let wtr = bcf::Writer::from_path("_test.vcf", &header, true, bcf::Format::Vcf).unwrap();
        unsafe { rust_htslib::htslib::bcf_hdr_sync(header.inner) };
        let header_t = unsafe { rust_htslib::htslib::bcf_hdr_dup(header.inner) };
        let mut header_view = HeaderView::new(header_t);

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
}
