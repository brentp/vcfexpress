use mlua::{AnyUserData, Lua, MetaMethod, UserData, UserDataFields, UserDataMethods, Value};
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
impl UserData for GTAllele {}
impl UserData for Genotype {}
impl UserData for Genotypes {}

use std::fmt;
impl fmt::Display for Genotype {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let Genotype(alleles) = self;
        write!(f, "{}", alleles[0].0)?;
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

pub fn register_genotypes(lua: &Lua) -> mlua::Result<()> {
    lua.register_userdata_type::<Genotype>(|reg| {
        reg.add_meta_function(MetaMethod::ToString, |_lua, this: AnyUserData| {
            let gts = format!("{}", this.borrow::<Genotype>()?);
            Ok(gts)
        });
        reg.add_field_method_get("alts", |_lua, this: &Genotype| {
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

        // index to get GTAllele
        reg.add_meta_function(
            MetaMethod::Index,
            |_lua, (this, idx): (AnyUserData, usize)| {
                let gts = this.borrow::<Genotype>()?;
                gts.0
                    .get(idx - 1)
                    .map(|allele| GTAllele(allele.0))
                    .ok_or_else(|| {
                        let msg =
                            format!("index out of bounds: {} in len: {}", idx - 1, gts.0.len());
                        mlua::Error::RuntimeError(msg)
                    })
            },
        );
    })?;

    lua.register_userdata_type::<GTAllele>(|reg| {
        reg.add_meta_function(MetaMethod::ToString, |_lua, this: AnyUserData| {
            Ok(this.borrow::<GTAllele>()?.0.to_string())
        });
        reg.add_field_method_get("phased", |_lua, this: &GTAllele| {
            Ok(match this.0 {
                GenotypeAllele::Phased(_) | GenotypeAllele::PhasedMissing => true,
                GenotypeAllele::Unphased(_) | GenotypeAllele::UnphasedMissing => false,
            })
        });
        reg.add_field_method_get("allele", |_lua, this: &GTAllele| {
            Ok(match this.0 {
                GenotypeAllele::Phased(i) => Value::Integer(i),
                GenotypeAllele::Unphased(i) => Value::Integer(i),
                _ => Value::Nil,
            })
        });
    })?;
    lua.register_userdata_type::<Genotypes>(|reg| {
        reg.add_meta_function(
            MetaMethod::Index,
            |_lua, (this, idx): (AnyUserData, usize)| {
                let ab = this.borrow::<Genotypes>()?;
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
        reg.add_meta_function(MetaMethod::Len, |_lua, this: AnyUserData| {
            let len = this.borrow::<Genotypes>()?.0.lock().0.len();
            Ok(len)
        });
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variant::{register_variant, Variant};
    use mlua::Lua;
    use rust_htslib::bcf;

    fn setup() -> (Lua, bcf::Record) {
        let lua = Lua::new();
        register_variant(&lua).expect("error registering variant");
        register_genotypes(&lua).expect("error registering genotypes");

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
        --print("printing from lua:", gts[i], "type:", type(i) )
        --print(gts[i][1], gts[i][2]) 
        --end
        local i = 1
        return tostring(gts[i])
        "#;
        let gts_exp = lua.load(gts_expr).set_name("gts").into_function().unwrap();
        let globals = lua.globals();
        let mut variant = Variant::new(record);

        lua.scope(|scope| {
            let ud = scope.create_any_userdata_ref_mut(&mut variant).unwrap();
            globals.raw_set("variant", ud).unwrap();
            let gtstring = gts_exp.call::<_, String>(());
            assert!(gtstring.is_ok());
            let gtstring = gtstring.unwrap();
            assert_eq!(gtstring, "0|1".to_string());

            // Add your assertions here...
            Ok(())
        })
        .unwrap();
    }
}
