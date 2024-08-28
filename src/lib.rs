//! This crate supports applying user-defined lua expressions to each variant in a VCF File.
//!
//pub mod genotypes;
//pub mod sample;
//pub mod header;
//pub mod pprint;
pub mod variant;
pub mod vcfexpress;

use rusty_v8 as v8;

pub fn register(isolate: &mut v8::Isolate, context: &v8::Local<v8::Context>) -> Result<(), Box<dyn std::error::Error>> {
    variant::register_variant(isolate, context)?;
    //header::register_header(isolate, context)?;
    //genotypes::register_genotypes(isolate, context)?;
    Ok(())
}
