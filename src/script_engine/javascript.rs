//! JavaScript engine implementation using rquickjs
//!
//! This module provides the JavaScript implementation of the ScriptEngine trait,
//! allowing users to write VCF filters and templates in JavaScript.

#[cfg(feature = "rquickjs")]
use super::{CompiledExpression, CompiledTemplate, ScriptConfig, ScriptEngine, ScriptError, ScriptResult};
#[cfg(feature = "rquickjs")]
use rust_htslib::bcf::header::HeaderView;
#[cfg(feature = "rquickjs")]
use rust_htslib::bcf::record::Record as VariantRecord;
#[cfg(feature = "rquickjs")]
use rquickjs::{Ctx, Runtime, Context, Object, Function, Result as JsResult, Value};
#[cfg(feature = "rquickjs")]
use std::collections::HashMap;
#[cfg(feature = "rquickjs")]
use crate::variant::Variant;

#[cfg(feature = "rquickjs")]
/// JavaScript implementation of the ScriptEngine trait
pub struct JSEngine {
    _runtime: Runtime,
    _ctx: Ctx<'static>,
    expressions: Vec<String>,
    templates: HashMap<String, String>,
    header: Option<HeaderView>,
}

#[cfg(feature = "rquickjs")]
impl ScriptEngine for JSEngine {
    fn new(config: &ScriptConfig) -> ScriptResult<Self> {
        let runtime = Runtime::new().map_err(ScriptError::JavaScript)?;
        let ctx = Ctx::full(&runtime).map_err(ScriptError::JavaScript)?;

        // Initialize JavaScript context
        ctx.with(|ctx| {
            // Load built-in utility functions
            let globals = ctx.globals();

            // Implement map, filter, all, any functions similar to Lua's
            ctx.eval::<(), _>(r#"
                // Utility functions for array/object operations
                globalThis.map = function(fn, arr, skipNil) {
                    if (Array.isArray(arr)) {
                        return arr.filter(x => skipNil && x == null ? false : true).map(fn);
                    } else {
                        const result = [];
                        for (const [key, value] of Object.entries(arr)) {
                            if (!skipNil || value != null) {
                                result.push(fn(value, key));
                            }
                        }
                        return result;
                    }
                };

                globalThis.filter = function(fn, arr, skipNil) {
                    if (Array.isArray(arr)) {
                        return arr.filter(x => skipNil && x == null ? false : true).filter(fn);
                    } else {
                        const result = [];
                        for (const [key, value] of Object.entries(arr)) {
                            if ((!skipNil || value != null) && fn(value, key)) {
                                result.push(value);
                            }
                        }
                        return result;
                    }
                };

                globalThis.all = function(fn, arr, skipNil) {
                    if (Array.isArray(arr)) {
                        return arr.filter(x => skipNil && x == null ? false : true).every(fn);
                    } else {
                        for (const value of Object.values(arr)) {
                            if (skipNil && value == null) continue;
                            if (!fn(value)) return false;
                        }
                        return true;
                    }
                };

                globalThis.any = function(fn, arr, skipNil) {
                    if (Array.isArray(arr)) {
                        return arr.filter(x => skipNil && x == null ? false : true).some(fn);
                    } else {
                        for (const value of Object.values(arr)) {
                            if (skipNil && value == null) continue;
                            if (fn(value)) return true;
                        }
                        return false;
                    }
                };

                // Pretty print function for debugging
                globalThis.pprint = function(obj, depth = 0) {
                    const indent = '  '.repeat(depth);
                    if (obj === null) return 'null';
                    if (obj === undefined) return 'undefined';
                    if (typeof obj !== 'object' || obj instanceof Date) {
                        return String(obj);
                    }
                    if (Array.isArray(obj)) {
                        if (obj.length === 0) return '[]';
                        return '[\n' + obj.map(v => indent + '  ' + pprint(v, depth + 1)).join(',\n') + '\n' + indent + ']';
                    }
                    const entries = Object.entries(obj);
                    if (entries.length === 0) return '{}';
                    return '{\n' + entries.map(([k, v]) =>
                        indent + '  ' + k + ' = ' + pprint(v, depth + 1)
                    ).join(',\n') + '\n' + indent + '}';
                };
            "#)?;

            // Apply sandbox if requested
            if config.sandbox {
                // In a real implementation, you would restrict access to dangerous APIs
                // For now, we'll just note that sandbox mode is enabled
                globals.set("SANDBOXED", true)?;
            }

            Ok::<(), rquickjs::Error>(())
        }).map_err(|e: rquickjs::Error| ScriptError::JavaScript(e))?;

        Ok(JSEngine {
            _runtime: runtime,
            _ctx: ctx,
            expressions: Vec::new(),
            templates: HashMap::new(),
            header: None,
        })
    }

    fn register_types(&mut self, header: &HeaderView) -> ScriptResult<()> {
        self.header = Some(header.clone());
        // In a real implementation, we would register Variant and Header types
        Ok(())
    }

    fn load_prelude(&mut self, prelude: &str) -> ScriptResult<()> {
        self._ctx.with(|ctx| {
            ctx.eval::<(), _>(prelude)
        }).map_err(|e: rquickjs::Error| ScriptError::JavaScript(e))?;
        Ok(())
    }

    fn compile_expression(&mut self, expr: &str) -> ScriptResult<CompiledExpression> {
        // Store as string for now
        self.expressions.push(expr.to_string());
        Ok(CompiledExpression::JavaScript(expr.to_string()))
    }

    fn compile_template(&mut self, template: &str) -> ScriptResult<CompiledTemplate> {
        // Store the template string for later evaluation
        self.templates.insert(template.to_string(), template.to_string());
        Ok(CompiledTemplate::JavaScript(template.to_string()))
    }

    fn evaluate_variant(&self, variant: &Variant, expr: &CompiledExpression) -> ScriptResult<bool> {
        match expr {
            CompiledExpression::JavaScript(expr_str) => {
                self._ctx.with(|ctx| {
                    // Create a proper JavaScript object that wraps the Variant record
                    let variant_obj = Self::create_variant_object(ctx, variant)?;

                    // Set the variant as global for the expression
                    let globals = ctx.globals();
                    globals.set("variant", variant_obj)?;

                    // Evaluate the expression
                    let expr_code = format!("return ({})", expr_str);
                    let result: bool = ctx.eval(expr_code.as_bytes())?;
                    Ok(result)
                }).map_err(|e: rquickjs::Error| ScriptError::JavaScript(e))
            }
            #[cfg(feature = "lua")]
            CompiledExpression::Lua(_) => Err(ScriptError::Message(
                "Cannot evaluate Lua expression with JavaScript engine".to_string()
            )),
        }
    }

    fn render_template(&self, variant: &Variant, template: &CompiledTemplate) -> ScriptResult<String> {
        match template {
            CompiledTemplate::JavaScript(template_str) => {
                self._ctx.with(|ctx| {
                    // Create a proper JavaScript object that wraps the Variant record
                    let variant_obj = Self::create_variant_object(ctx, variant)?;

                    // Set the variant as global for the template
                    let globals = ctx.globals();
                    globals.set("variant", variant_obj)?;

                    // For JavaScript templates, we expect template literals
                    // If it's not already a template literal, wrap it in backticks
                    let eval_template = if template_str.contains('`') {
                        template_str.clone()
                    } else {
                        format!("`{}`", template_str)
                    };

                    // Evaluate the template
                    let result: String = ctx.eval(eval_template.as_bytes())?;
                    Ok(result)
                }).map_err(|e: rquickjs::Error| ScriptError::JavaScript(e))
            }
            #[cfg(feature = "lua")]
            CompiledTemplate::Lua(_) => Err(ScriptError::Message(
                "Cannot render Lua template with JavaScript engine".to_string()
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
            CompiledExpression::JavaScript(expr_str) => {
                self._ctx.with(|ctx| {
                    // Create simple variant context in JavaScript
                    let pos = variant.pos();
                    let qual = variant.qual();

                    // Create variant object in JavaScript
                    let js_code = format!(r#"
                        globalThis.variant = {{
                            pos: {},
                            qual: {}
                        }};
                    "#, pos, qual);
                    ctx.eval::<(), _>(js_code.as_bytes())?;

                    // Evaluate the expression as a string and get the result
                    let expr_code = format!("String({})", expr_str);
                    let result_str: String = ctx.eval(expr_code.as_bytes())?;

                    // Try to parse the result as different types
                    // For now, we'll handle basic cases - this could be enhanced
                    if result_str == "true" || result_str == "false" {
                        // Boolean flag
                        let flag_value = result_str == "true";
                        let field_bytes = field.as_bytes();
                        if flag_value {
                            variant.push_info_flag(field_bytes)
                                .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, format!("HTSlib error: {}", e)))?;
                        } else {
                            variant.clear_info_flag(field_bytes)
                                .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, format!("HTSlib error: {}", e)))?;
                        }
                    } else if let Ok(int_val) = result_str.parse::<i32>() {
                        // Integer value
                        variant.push_info_integer(field.as_bytes(), &[int_val])
                            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, format!("HTSlib error: {}", e)))?;
                    } else if let Ok(float_val) = result_str.parse::<f32>() {
                        // Float value
                        variant.push_info_float(field.as_bytes(), &[float_val])
                            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, format!("HTSlib error: {}", e)))?;
                    } else {
                        // String value (skip empty strings)
                        if !result_str.is_empty() && result_str != "null" && result_str != "undefined" {
                            variant.push_info_string(field.as_bytes(), &[result_str.as_bytes()])
                                .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, format!("HTSlib error: {}", e)))?;
                        }
                    }

                    Ok(())
                }).map_err(|e: rquickjs::Error| ScriptError::JavaScript(e))
            }
            #[cfg(feature = "lua")]
            CompiledExpression::Lua(_) => Err(ScriptError::Message(
                "Cannot evaluate Lua expression with JavaScript engine".to_string()
            )),
        }
    }

    fn language(&self) -> super::ScriptLanguage {
        super::ScriptLanguage::JavaScript
    }

    /// Create a JavaScript object that wraps the Variant record using JavaScript code injection
    /// This is a more practical approach given rquickjs API limitations
    fn create_variant_object_js(variant: &Variant) -> String {
        let record = variant.record();
        let header = record.header();

        // Get basic variant information
        let chrom = record.rid()
            .and_then(|id| header.rid2name(id).ok())
            .and_then(|name| std::str::from_utf8(name).ok())
            .unwrap_or("unknown");
        let pos = record.pos();
        let end = record.end();
        let id_bytes = record.id();
        let id = std::str::from_utf8(&id_bytes).unwrap_or(".");
        let qual = record.qual();

        // Get alleles
        let alleles = record.alleles();
        let ref_allele = std::str::from_utf8(alleles[0]).unwrap_or("");
        let alt_alleles: Vec<String> = alleles.iter()
            .skip(1)
            .map(|&allele| std::str::from_utf8(allele).unwrap_or("").to_string())
            .collect();

        // Get filters
        let filter_names: Vec<String> = record.filters()
            .iter()
            .map(|&id| std::str::from_utf8(&header.id_to_name(id)).unwrap_or("").to_string())
            .collect();

        // Create JavaScript code that builds the variant object
        format!(r#"
        globalThis.variant = {{
            chrom: "{chrom}",
            pos: {pos},
            start: {pos},
            stop: {end},
            id: "{id}",
            qual: {qual},
            REF: "{ref_allele}",
            ALT: [{alt_js}],
            filters: [{filter_js}],
            ref_allele: "{ref_allele}",
            alt_alleles: [{alt_vec_js}],
            alleles: [{alleles_js_vec}],
            info: function(key) {{
                // Simplified info function - returns undefined for now
                // In a full implementation, this would access the actual INFO fields
                return undefined;
            }},
            format: function(key) {{
                // Simplified format function - returns undefined for now
                // In a full implementation, this would access the actual FORMAT fields
                return undefined;
            }}
        }};
        "#,
            chrom = chrom,
            pos = pos,
            end = end,
            id = id,
            qual = qual,
            ref_allele = ref_allele,
            alt_js = alt_alleles.iter().map(|alt| format!("'{}'", alt)).collect::<Vec<_>>().join(", "),
            filter_js = filter_names.iter().map(|f| format!("'{}'", f)).collect::<Vec<_>>().join(", "),
            alt_vec_js = alt_alleles.iter().map(|alt| format!("'{}'", alt)).collect::<Vec<_>>().join(", "),
            alleles_js_vec = alleles.iter()
                .map(|&allele| format!("'{}'", std::str::from_utf8(allele).unwrap_or("")))
                .collect::<Vec<_>>().join(", ")
        )
    }
}