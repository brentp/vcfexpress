pub mod genotypes;
//pub mod sample;
pub mod variant;

pub fn register(lua: &mlua::Lua) -> mlua::Result<()> {
    variant::register_variant(lua)?;
    genotypes::register_genotypes(lua)?;
    Ok(())
}
