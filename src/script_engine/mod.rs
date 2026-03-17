//! Abstraction layer for supporting multiple scripting languages
//!
//! This module defines the common interface that all scripting engines must implement,
//! allowing VCFExpress to support different languages like Lua and JavaScript.

use rust_htslib::bcf::header::HeaderView;
use rust_htslib::bcf::record::Record as VariantRecord;
use std::fmt;
use crate::variant::Variant;

/// Configuration for scripting engines
#[derive(Debug, Clone)]
pub struct ScriptConfig {
    /// Whether to run scripts in a sandboxed environment
    pub sandbox: bool,
    /// The specific language to use
    pub language: ScriptLanguage,
}

/// Supported scripting languages
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScriptLanguage {
    /// Lua/Luau scripting (default, maintains backward compatibility)
    Lua,
    /// JavaScript using rquickjs
    JavaScript,
}

impl Default for ScriptLanguage {
    fn default() -> Self {
        ScriptLanguage::Lua
    }
}

impl fmt::Display for ScriptLanguage {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ScriptLanguage::Lua => write!(f, "lua"),
            ScriptLanguage::JavaScript => write!(f, "javascript"),
        }
    }
}

impl std::str::FromStr for ScriptLanguage {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "lua" => Ok(ScriptLanguage::Lua),
            "javascript" | "js" => Ok(ScriptLanguage::JavaScript),
            _ => Err(format!("Unsupported scripting language: {}", s)),
        }
    }
}

/// Compiled expression ready for evaluation
#[derive(Debug)]
pub enum CompiledExpression {
    #[cfg(feature = "lua")]
    Lua(mlua::Function),
    #[cfg(feature = "javascript")]
    JavaScript(String), // Store as string for now, will be compiled when needed
}

/// Compiled template ready for rendering
#[derive(Debug)]
pub enum CompiledTemplate {
    #[cfg(feature = "lua")]
    Lua(mlua::Function),
    #[cfg(feature = "javascript")]
    JavaScript(String), // Store as string for now, will be compiled when needed
}

/// Errors that can occur during scripting operations
#[derive(Debug)]
pub enum ScriptError {
    /// Error from Lua engine
    #[cfg(feature = "lua")]
    Lua(mlua::Error),
    /// Error from JavaScript engine
    #[cfg(feature = "javascript")]
    JavaScript(rquickjs::Error),
    /// Language not supported (feature not enabled)
    LanguageNotSupported(ScriptLanguage),
    /// Generic error message
    Message(String),
}

impl fmt::Display for ScriptError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            #[cfg(feature = "lua")]
            ScriptError::Lua(e) => write!(f, "Lua error: {}", e),
            #[cfg(feature = "javascript")]
            ScriptError::JavaScript(e) => write!(f, "JavaScript error: {}", e),
            ScriptError::LanguageNotSupported(lang) => {
                write!(f, "Language '{}' not supported in this build", lang)
            }
            ScriptError::Message(msg) => write!(f, "Script error: {}", msg),
        }
    }
}

impl std::error::Error for ScriptError {}

/// Result type for scripting operations
pub type ScriptResult<T> = Result<T, ScriptError>;

/// Common interface for all scripting engines
///
/// This trait allows VCFExpress to support multiple scripting languages
/// while keeping the core logic language-agnostic.
pub trait ScriptEngine {
    /// Create a new scripting engine with the given configuration
    fn new(config: &ScriptConfig) -> ScriptResult<Self>
    where
        Self: Sized;

    /// Register types (Variant, Header, etc.) with the scripting engine
    /// This is called once before processing begins
    fn register_types(&mut self, header: &HeaderView) -> ScriptResult<()>;

    /// Load and execute prelude code (runs once before any variants are processed)
    fn load_prelude(&mut self, prelude: &str) -> ScriptResult<()>;

    /// Compile a boolean expression for filtering
    fn compile_expression(&mut self, expr: &str) -> ScriptResult<CompiledExpression>;

    /// Compile a template string for output formatting
    fn compile_template(&mut self, template: &str) -> ScriptResult<CompiledTemplate>;

    /// Evaluate a boolean expression against a variant
    fn evaluate_variant(&self, variant: &Variant, expr: &CompiledExpression) -> ScriptResult<bool>;

    /// Render a template with the variant's data
    fn render_template(&self, variant: &Variant, template: &CompiledTemplate) -> ScriptResult<String>;

    /// Set an INFO field on a variant record based on an expression
    fn set_info_field(
        &self,
        variant: &mut VariantRecord,
        field: &str,
        expr: &CompiledExpression,
    ) -> ScriptResult<()>;

    /// Get the language this engine supports
    fn language(&self) -> ScriptLanguage;
}

/// Boxed trait object for scripting engines
pub type BoxedEngine = Box<dyn ScriptEngine>;

/// Factory function to create a script engine based on configuration
pub fn create_engine(config: &ScriptConfig) -> ScriptResult<BoxedEngine> {
    match config.language {
        #[cfg(feature = "lua")]
        ScriptLanguage::Lua => {
            let engine = lua::LuaEngine::new(config)?;
            Ok(Box::new(engine))
        }
        #[cfg(feature = "javascript")]
        ScriptLanguage::JavaScript => {
            let engine = javascript::JSEngine::new(config)?;
            Ok(Box::new(engine))
        }
        #[cfg(not(feature = "lua"))]
        ScriptLanguage::Lua => {
            Err(ScriptError::LanguageNotSupported(ScriptLanguage::Lua))
        }
        #[cfg(not(feature = "javascript"))]
        ScriptLanguage::JavaScript => {
            Err(ScriptError::LanguageNotSupported(ScriptLanguage::JavaScript))
        }
    }
}

#[cfg(feature = "lua")]
pub mod lua;

#[cfg(feature = "javascript")]
pub mod javascript;