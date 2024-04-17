pub mod genotypes;
//pub mod sample;
pub mod variant;

pub fn register(lua: &mlua::Lua) -> mlua::Result<()> {
    genotypes::register_genotypes(lua)?;
    variant::register_variant(lua)?;
    Ok(())
}
