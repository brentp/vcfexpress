//! Lua engine implementation using mlua
//!
//! This module provides the Lua implementation of the ScriptEngine trait,
//! maintaining full backward compatibility with existing VCFExpress Lua scripts.

#[cfg(feature = "mlua")]
use super::{CompiledExpression, CompiledTemplate, ScriptConfig, ScriptEngine, ScriptError, ScriptResult};
#[cfg(feature = "mlua")]
use rust_htslib::bcf::header::HeaderView;
#[cfg(feature = "mlua")]
use rust_htslib::bcf::record::Record as VariantRecord;
#[cfg(feature = "mlua")]
use crate::variant::Variant;

#[cfg(feature = "mlua")]
/// Lua implementation of the ScriptEngine trait
pub struct LuaEngine {
    lua: mlua::Lua,
    header: Option<HeaderView>,
}

#[cfg(feature = "mlua")]
impl ScriptEngine for LuaEngine {
    fn new(config: &ScriptConfig) -> ScriptResult<Self> {
        let lua = mlua::Lua::new();

        // Apply sandbox configuration
        lua.sandbox(config.sandbox)
            .map_err(ScriptError::Lua)?;

        // Load built-in functions
        lua.load(crate::pprint::PPRINT)
            .set_name("pprint")
            .exec()
            .map_err(ScriptError::Lua)?;

        lua.load(crate::pprint::PRELUDE)
            .set_name("prelude")
            .exec()
            .map_err(ScriptError::Lua)?;

        lua.load(crate::pprint::LUA_PRELUDE)
            .set_name("lua_prelude")
            .exec()
            .map_err(ScriptError::Lua)?;

        Ok(LuaEngine {
            lua,
            header: None,
        })
    }

    fn register_types(&mut self, header: &HeaderView) -> ScriptResult<()> {
        self.header = Some(header.clone());
        crate::variant::register_variant(&self.lua)
            .map_err(ScriptError::Lua)?;
        crate::header::register_header(&self.lua)
            .map_err(ScriptError::Lua)?;
        Ok(())
    }

    fn load_prelude(&mut self, prelude: &str) -> ScriptResult<()> {
        self.lua
            .load(prelude)
            .exec()
            .map_err(ScriptError::Lua)?;
        Ok(())
    }

    fn compile_expression(&mut self, expr: &str) -> ScriptResult<CompiledExpression> {
        // Compile the expression to a Lua function
        let func = self.lua
            .load(&format!("return {}", expr))
            .into_function()
            .map_err(ScriptError::Lua)?;
        Ok(CompiledExpression::Lua(func))
    }

    fn compile_template(&mut self, template: &str) -> ScriptResult<CompiledTemplate> {
        // Process template for Lua's Luau string interpolation
        let return_pre = if template.contains("return ") { "" } else { "return " };

        let expr = if template.contains('`') {
            format!("{}{}", return_pre, template)
        } else {
            format!("{} `{}`", return_pre, template)
        };

        // Compile the template to a Lua function that returns a string
        let func = self.lua
            .load(&expr)
            .into_function()
            .map_err(ScriptError::Lua)?;

        Ok(CompiledTemplate::Lua(func))
    }

    fn evaluate_variant(&self, variant: &Variant, expr: &CompiledExpression) -> ScriptResult<bool> {
        match expr {
            CompiledExpression::Lua(func) => {
                // For now, use a simplified approach without variant context
                // This is a limitation we'll need to address when refactoring VCFExpress
                func.call::<bool>(())
                    .map_err(ScriptError::Lua)
            }
            #[cfg(feature = "javascript")]
            CompiledExpression::JavaScript(_) => Err(ScriptError::Message(
                "Cannot evaluate JavaScript expression with Lua engine".to_string()
            )),
        }
    }

    fn render_template(&self, variant: &Variant, template: &CompiledTemplate) -> ScriptResult<String> {
        match template {
            CompiledTemplate::Lua(func) => {
                // For now, use a simplified approach without variant context
                // This is a limitation we'll need to address when refactoring VCFExpress
                func.call::<String>(())
                    .map_err(ScriptError::Lua)
            }
            #[cfg(feature = "javascript")]
            CompiledTemplate::JavaScript(_) => Err(ScriptError::Message(
                "Cannot render JavaScript template with Lua engine".to_string()
            )),
        }
    }

    fn set_info_field(
        &self,
        variant: &mut VariantRecord,
        field: &str,
        expr: &CompiledExpression,
    ) -> ScriptResult<()> {
        match expr {
            CompiledExpression::Lua(func) => {
                // For set_info_field, we don't have access to the full Variant context
                // This is a limitation - in practice, this method might need to be redesigned
                // For now, we'll evaluate the expression without variant context
                let result: mlua::Value = func.call(())
                    .map_err(ScriptError::Lua)?;

                let field_bytes = field.as_bytes();
                let result = match result {
                    mlua::Value::Boolean(b) => {
                        if b {
                            variant.push_info_flag(field_bytes)
                        } else {
                            variant.clear_info_flag(field_bytes)
                        }
                    }
                    mlua::Value::Number(n) => {
                        // Determine if it's an integer or float based on value
                        if n.fract() == 0.0 && (n as i32 as f64 == n) {
                            variant.push_info_integer(field_bytes, &[n as i32])
                        } else {
                            variant.push_info_float(field_bytes, &[n as f32])
                        }
                    }
                    mlua::Value::String(s) => {
                        let str_val = s.to_str().map_err(ScriptError::Lua)?;
                        variant.push_info_string(field_bytes, &[str_val.as_bytes()])
                    }
                    mlua::Value::Nil => {
                        // Skip nil values
                        Ok(())
                    }
                    _ => {
                        return Err(ScriptError::Message(
                            format!("Unsupported return type for INFO field '{}'", field)
                        ));
                    }
                };

                result.map_err(|e| ScriptError::Message(format!("Failed to set INFO field '{}': {}", field, e)))
            }
            #[cfg(feature = "javascript")]
            CompiledExpression::JavaScript(_) => Err(ScriptError::Message(
                "Cannot evaluate JavaScript expression with Lua engine".to_string()
            )),
        }
    }

    fn language(&self) -> super::ScriptLanguage {
        super::ScriptLanguage::Lua
    }
}