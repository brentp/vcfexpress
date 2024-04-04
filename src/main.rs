use fasteval::Evaler;
use crate::VariantRecord;
use std::collections::HashMap;

/// Filters `VariantRecord` instances based on FastEval expressions.
pub fn filter_records<'a, I>(
    records: I,
    expressions: &[&str],
) -> impl Iterator<Item = VariantRecord> + 'a
where
    I: Iterator<Item = VariantRecord> + 'a,
{
    records.filter(move |record| {
        expressions.iter().all(|&expr| {
            let mut slab = fasteval::Slab::new();
            let mut ns = HashMap::new();
            for (key, value) in &record.info {
                ns.insert(format!("INFO__{}", key), value.parse::<f64>().unwrap_or(0.0));
            }
            fasteval::ez_eval(expr, &mut slab, &mut ns).unwrap_or(0.0) != 0.0
        })
    })
}