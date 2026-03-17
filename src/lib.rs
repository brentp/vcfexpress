//! This crate supports applying user-defined lua or javascript expressions to each variant in a VCF File.
//!
pub mod genotypes;
//pub mod sample;
pub mod header;
pub mod pprint;
pub mod variant;
pub mod vcfexpress;
pub mod script_engine;
pub mod error;

// Legacy register function for backward compatibility
#[cfg(feature = "lua")]
pub fn register(lua: &mlua::Lua) -> mlua::Result<()> {
    variant::register_variant(lua)?;
    header::register_header(lua)
}
