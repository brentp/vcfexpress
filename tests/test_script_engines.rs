//! Tests for the ScriptEngine abstraction layer and language implementations

#[cfg(test)]
mod tests {
    use vcfexpress::script_engine::{ScriptConfig, ScriptLanguage, create_engine};

    #[test]
    #[cfg(feature = "lua")]
    fn test_lua_engine_creation() {
        let config = ScriptConfig {
            sandbox: false,
            language: ScriptLanguage::Lua,
        };

        let engine = create_engine(&config);
        assert!(engine.is_ok(), "Failed to create Lua engine");
    }

    #[test]
    #[cfg(feature = "javascript")]
    fn test_javascript_engine_creation() {
        let config = ScriptConfig {
            sandbox: false,
            language: ScriptLanguage::JavaScript,
        };

        let engine = create_engine(&config);
        assert!(engine.is_ok(), "Failed to create JavaScript engine");
    }

    #[test]
    #[cfg(not(feature = "lua"))]
    fn test_lua_engine_disabled() {
        let config = ScriptConfig {
            sandbox: false,
            language: ScriptLanguage::Lua,
        };

        let engine = create_engine(&config);
        assert!(engine.is_err(), "Lua engine should not be available without lua feature");
    }

    #[test]
    #[cfg(not(feature = "javascript"))]
    fn test_javascript_engine_disabled() {
        let config = ScriptConfig {
            sandbox: false,
            language: ScriptLanguage::JavaScript,
        };

        let engine = create_engine(&config);
        assert!(engine.is_err(), "JavaScript engine should not be available without javascript feature");
    }

    #[test]
    fn test_script_language_parsing() {
        assert_eq!("lua".parse::<ScriptLanguage>().unwrap(), ScriptLanguage::Lua);
        assert_eq!("lua".parse::<ScriptLanguage>().unwrap(), "LUA".parse::<ScriptLanguage>().unwrap());
        assert_eq!("javascript".parse::<ScriptLanguage>().unwrap(), ScriptLanguage::JavaScript);
        assert_eq!("js".parse::<ScriptLanguage>().unwrap(), ScriptLanguage::JavaScript);
        assert!("python".parse::<ScriptLanguage>().is_err());
    }

    #[test]
    fn test_script_language_display() {
        assert_eq!(format!("{}", ScriptLanguage::Lua), "lua");
        assert_eq!(format!("{}", ScriptLanguage::JavaScript), "javascript");
    }
}