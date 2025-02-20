use mlua::{MetaMethod, UserData, UserDataFields, UserDataMethods, Value};
use parking_lot::Mutex;
use rust_htslib::bcf;
use rust_htslib::bcf::record::{self, GenotypeAllele};
use std::sync::Arc;

pub(crate) struct I32Buffer(
    pub(crate) bcf::record::BufferBacked<'static, Vec<&'static [i32]>, record::Buffer>,
);

struct GTAllele(bcf::record::GenotypeAllele);
struct Genotype(Vec<GTAllele>);

pub(crate) struct Genotypes(pub(crate) Arc<Mutex<I32Buffer>>);

impl std::fmt::Debug for GTAllele {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self.0)
    }
}

unsafe impl Send for I32Buffer {}
impl UserData for I32Buffer {}
impl UserData for GTAllele {
    fn add_methods<M: UserDataMethods<Self>>(methods: &mut M) {
        methods.add_meta_method(MetaMethod::ToString, |_lua, this, ()| {
            Ok(this.0.to_string())
        });
    }

    fn add_fields<M: UserDataFields<Self>>(fields: &mut M) {
        fields.add_field_method_get("phased", |_lua, this| {
            Ok(match this.0 {
                GenotypeAllele::Phased(_) | GenotypeAllele::PhasedMissing => true,
                GenotypeAllele::Unphased(_) | GenotypeAllele::UnphasedMissing => false,
            })
        });
        fields.add_field_method_get("allele", |_lua, this| {
            Ok(match this.0 {
                GenotypeAllele::Phased(i) => Value::Integer(i),
                GenotypeAllele::Unphased(i) => Value::Integer(i),
                _ => Value::Nil,
            })
        });
    }
}

impl UserData for Genotype {
    fn add_fields<M: UserDataFields<Self>>(fields: &mut M) {
        fields.add_field_method_get("alts", |_lua, this| {
            Ok(this
                .0
                .iter()
                .map(|x| match x.0 {
                    GenotypeAllele::Phased(i) => i,
                    GenotypeAllele::Unphased(i) => i,
                    GenotypeAllele::PhasedMissing => 0,
                    GenotypeAllele::UnphasedMissing => 0,
                })
                .sum::<i32>())
        })
    }
    fn add_methods<M: UserDataMethods<Self>>(methods: &mut M) {
        methods.add_meta_method(MetaMethod::ToString, |_lua, this, ()| {
            let gts = format!("{}", this);
            Ok(gts)
        });

        methods.add_method("alts", |_lua, this, ()| {
            Ok(this
                .0
                .iter()
                .map(|x| match x.0 {
                    GenotypeAllele::Phased(i) => {
                        if i != 0 {
                            1
                        } else {
                            0
                        }
                    }
                    GenotypeAllele::Unphased(i) => {
                        if i != 0 {
                            1
                        } else {
                            0
                        }
                    }
                    _ => 0,
                })
                .sum::<i32>())
        });

        methods.add_meta_method(MetaMethod::Index, |_lua, this, idx: usize| {
            this.0
                .get(idx - 1)
                .map(|allele| GTAllele(allele.0))
                .ok_or_else(|| {
                    let msg = format!("index out of bounds: {} in len: {}", idx - 1, this.0.len());
                    mlua::Error::RuntimeError(msg)
                })
        });
    }
}
impl UserData for Genotypes {
    fn add_methods<M: UserDataMethods<Self>>(methods: &mut M) {
        methods.add_meta_method(
            MetaMethod::Index,
            |_lua, this, idx: usize| -> mlua::Result<Genotype> {
                let ab = this;
                let buffer = &ab.0.lock().0;
                let len = buffer.len();
                buffer
                    .get(idx - 1)
                    .map(|&x| {
                        let gts = x
                            .iter()
                            .map(|&allele_int| {
                                GTAllele(bcf::record::GenotypeAllele::from(allele_int))
                            })
                            .collect::<Vec<GTAllele>>();
                        Genotype(gts)
                    })
                    .ok_or_else(|| {
                        let msg = format!("index out of bounds: {} in len: {}", idx - 1, len);
                        mlua::Error::RuntimeError(msg)
                    })
            },
        );

        methods.add_meta_method(MetaMethod::Len, |_lua, this, ()| -> mlua::Result<usize> {
            let len = this.0.lock().0.len();
            Ok(len)
        });
    }
}

use std::fmt;
impl fmt::Display for Genotype {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let Genotype(alleles) = self;
        write!(f, "{}", alleles[0].0)?;
        // convert to the alleles
        for allele in alleles[1..].iter() {
            let allele = allele.0;
            let sep = match allele {
                GenotypeAllele::Phased(_) | GenotypeAllele::PhasedMissing => "|",
                GenotypeAllele::Unphased(_) | GenotypeAllele::UnphasedMissing => "/",
            };
            write!(f, "{}{}", sep, allele)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::variant::{register_variant, HeaderMap, Variant};
    use mlua::Lua;
    use rust_htslib::bcf;

    fn setup() -> (Lua, bcf::Record) {
        let lua = Lua::new();
        register_variant(&lua).expect("error registering variant");
        let mut header = bcf::Header::new();
        header.push_record(r#"##contig=<ID=chr1,length=10000>"#.as_bytes());
        header.push_record(
            r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#.as_bytes(),
        );
        header.push_sample("NA12878".as_bytes());
        header.push_sample("NA12879".as_bytes());
        let tmp_path = "_test.bcf";
        let vcf = bcf::Writer::from_path(tmp_path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let _ = record.set_rid(Some(vcf.header().name2rid(b"chr1").unwrap()));
        record.set_pos(6);
        record.set_id(b"rs1234").unwrap();
        let alleles = &[
            bcf::record::GenotypeAllele::Unphased(0),
            bcf::record::GenotypeAllele::Phased(1),
            bcf::record::GenotypeAllele::Unphased(1),
            bcf::record::GenotypeAllele::Unphased(1),
        ];
        record.push_genotypes(alleles).unwrap();

        (lua, record)
    }

    #[test]
    fn test_gts_expression() {
        let (lua, record) = setup();
        let gts_expr = r#"local gts = variant.genotypes; 
        --for i = 1, #gts do 
        --  print("printing from lua:", gts[i], "type:", type(i) )
        -- print(gts[i][1], gts[i][2]) 
        --end
        local i = 1
        local gt = gts[i]
        return tostring(gt) .. " " .. tostring(gt[1]) .. " " .. tostring(gt[2]) .. " " .. tostring(gt.alts)  .. " " .. tostring(gt[2].phased) .. " "  .. tostring(gt[2].allele)
        "#;
        let gts_exp = lua.load(gts_expr).set_name("gts").into_function().unwrap();
        let globals = lua.globals();
        let mut variant = Variant::new(record, HeaderMap::new());

        lua.scope(|scope| {
            let ud = scope.create_any_userdata_ref_mut(&mut variant).unwrap();
            globals.raw_set("variant", ud).unwrap();
            let gtstring = gts_exp.call::<String>(());
            eprintln!("gtstring: {:?}", gtstring);
            assert!(gtstring.is_ok());
            let gtstring = gtstring.unwrap();
            assert_eq!(gtstring, "0|1 0 1 1 true 1".to_string());

            // Add your assertions here...
            Ok(())
        })
        .unwrap();
    }
}
