//! Unified error handling for VCFExpress
//!
//! This module consolidates error types from different components
//! including scripting engines, VCF processing, and I/O operations.

use crate::script_engine::ScriptError;
use std::fmt;

/// Main error type for VCFExpress
#[derive(Debug)]
pub enum VCFError {
    /// Error during VCF/BCF file operations
    Htslib(rust_htslib::errors::Error),
    /// Error during scripting operations
    Script(ScriptError),
    /// I/O error
    Io(std::io::Error),
    /// Error with command line arguments
    Cli(String),
    /// Error during template rendering
    Template(String),
    /// Error during info field operations
    InfoField(String),
}

impl fmt::Display for VCFError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            VCFError::Htslib(e) => write!(f, "HTSlib error: {}", e),
            VCFError::Script(e) => write!(f, "Script error: {}", e),
            VCFError::Io(e) => write!(f, "I/O error: {}", e),
            VCFError::Cli(msg) => write!(f, "CLI error: {}", msg),
            VCFError::Template(msg) => write!(f, "Template error: {}", msg),
            VCFError::InfoField(msg) => write!(f, "Info field error: {}", msg),
        }
    }
}

impl std::error::Error for VCFError {}

impl From<rust_htslib::errors::Error> for VCFError {
    fn from(err: rust_htslib::errors::Error) -> Self {
        VCFError::Htslib(err)
    }
}

impl From<std::io::Error> for VCFError {
    fn from(err: std::io::Error) -> Self {
        VCFError::Io(err)
    }
}

impl From<ScriptError> for VCFError {
    fn from(err: ScriptError) -> Self {
        VCFError::Script(err)
    }
}

/// Result type for VCFExpress operations
pub type VCFResult<T> = Result<T, VCFError>;