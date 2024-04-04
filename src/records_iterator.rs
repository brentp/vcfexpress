use std::collections::HashMap;

/// Represents a single variant record from a VCF/BCF file.
pub struct VariantRecord {
    /// Stores values from the INFO field.
    pub info: HashMap<String, String>,
}
