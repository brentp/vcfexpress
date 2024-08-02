//! This crate supports applying user-defined lua expressions to each variant in a VCF File.
//!
pub mod genotypes;
//pub mod sample;
pub mod header;
pub mod pprint;
pub mod variant;
pub mod vcfexpress;

pub fn register(lua: &mlua::Lua) -> mlua::Result<()> {
    variant::register_variant(lua)?;
    genotypes::register_genotypes(lua)?;
    header::register_header(lua)
}
