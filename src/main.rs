use fasteval::Evaler;
use crate::VariantRecord;
use std::collections::HashMap;

use mlua::{Lua, Value};
use crate::VariantRecord;
use std::collections::HashMap;

/// Evaluates Lua expressions on a VariantRecord's INFO field and returns whether it matches the filter.
fn evaluate_expression(lua: &Lua, record: &VariantRecord, expression: &str) -> Result<bool, mlua::Error> {
    let lua_exp = format!("return {}", expression);
    let func: mlua::Function = lua.load(&lua_exp).eval()?;
    
    // Convert VariantRecord INFO field to Lua table
    let info_table = lua.create_table()?;
    for (key, value) in &record.info {
        info_table.set(key, value)?;
    }

    func.call::<_, bool>((info_table,))
}

/// Filters variants based on user-provided Lua expressions.
pub fn filter_variants<'a, I>(
    variants: I,
    expressions: &[&str],
) -> Result<impl Iterator<Item = VariantRecord> + 'a, mlua::Error>
where
    I: Iterator<Item = VariantRecord> + 'a,
{
    let lua = Lua::new();
    Ok(variants.filter(move |variant| {
        expressions.iter().all(|&expression| {
            evaluate_expression(&lua, variant, expression).unwrap_or(false)
        })
    }))
}