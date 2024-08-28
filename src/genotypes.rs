use rusty_v8 as v8;
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

pub fn register_genotypes(isolate: &mut v8::Isolate, context: &v8::Local<v8::Context>) -> Result<(), Box<dyn std::error::Error>> {
    let scope = &mut v8::HandleScope::with_context(isolate, context);
    let global = context.global(scope);

    let genotype_template = v8::ObjectTemplate::new(scope);
    
    // Add methods
    genotype_template.set(
        v8::String::new(scope, "toString").unwrap().into(),
        v8::FunctionTemplate::new(scope, genotype_to_string).into(),
    );


    // expose so we can get new Genotype() objects
    global.set(
        scope,
        v8::String::new(scope, "Genotype").unwrap().into(),
        genotype_template.into(),
    )?
}

// Implement methods
fn genotype_to_string(
    scope: &mut v8::HandleScope,
    args: v8::FunctionCallbackArguments,
    mut retval: v8::ReturnValue,
) {
    let genotype = args.this().get_internal_field(scope, 0).unwrap().is_external().unwrap();
    let genotype = unsafe { &*(genotype.value() as *const Genotype) };
    
    let gts = format!("{}", genotype);
    retval.set(v8::String::new(scope, &gts).unwrap().into());
}

// ... implement other methods

// Update the tests to use V8 instead of Lua
#[cfg(test)]
mod tests {
    use super::*;
    use v8;

    fn setup() -> (v8::Isolate, bcf::Record) {
        let isolate = v8::Isolate::new(v8::CreateParams::default());
        let handle_scope = &mut v8::HandleScope::new(&mut isolate);
        let context = v8::Context::new(handle_scope);

        register_genotypes(&mut isolate, &context).expect("error registering genotypes");

        // Create a dummy bcf::Record
        let header = bcf::Header::new();
        let mut record = bcf::Record::new(&header);

        (isolate, record)
    }

    #[test]
    fn test_gts_expression() {
        let (isolate, record) = setup();
        let scope = &mut v8::HandleScope::new(&isolate);
        let context = v8::Context::new(scope);

        let gts_expr = r#"
        let gts = variant.genotypes; 
        return gts[1].toString();
        "#;

        let script = v8::Script::compile(scope, v8::String::new(scope, gts_expr).unwrap(), None).unwrap();
        let result = script.run(scope).unwrap();
        let gtstring = result.to_string(scope).unwrap().to_rust_string_lossy(scope);

        assert_eq!(gtstring, "0|1");
    }
}
