# VCFExpress Lua API Documentation

This document details all Lua attributes and functions available when using VCFExpress.

## Global Functions

These functions are available globally for array/table operations:

- `map(function, table, [skip_nil])`: Apply a function to each element in a table
- `filter(function, table, [skip_nil])`: Filter table elements based on a predicate function
- `all(function, table, [skip_nil])`: Check if all elements satisfy a predicate function
- `any(function, table, [skip_nil])`: Check if any element satisfies a predicate function
- `pprint(table)`: Pretty print a table structure

- `skip_nil` in the functions above means that values in the table that are `nil` will be skipped.
- *NOTE* that for iteration, all of these functions will use `ipairs` if `#tbl > 0` and `pairs` otherwise.

## Variant Object

The main variant object provides access to VCF record data:

### Basic Attributes

- `variant.chrom` (string): Chromosome name
- `variant.pos` (integer, get/set): 0-based position
- `variant.start` (integer): Start position (same as pos)
- `variant.stop` (integer): End position
- `variant.qual` (number, get/set): Variant quality score
- `variant.id` (string, get/set): Variant ID
- `variant.REF` (string, get/set): Reference allele
- `variant.ALT` (table, get/set): Array of alternate alleles
- `variant.FILTER` (string, get/set): First filter value
- `variant.filters` (table, get/set): Array of all filter values
- `variant.genotypes`: Array of genotype objects for all samples

### Methods

- `variant:info(field_name, [index])`: Get INFO field value(s)
  - Returns number/string/bool/table depending on field type
  - Optional index parameter to get specific value for multi-value fields
- `variant:format(field_name)`: Get FORMAT field values for all samples
  - Returns array of values, type depends on field definition
- `variant:sample(sample_name)`: Get all FORMAT fields for a specific sample
  - Returns table with format field values
  - Special handling for GT field to provide allele and phase information

## Genotype Objects

Accessed through `variant.genotypes[sample_index]`:

### Attributes

- `genotype.alts`: Number of alternate alleles (non-reference, non-missing)
- `tostring(genotype)`: String representation (e.g. "0/1" or "1|0")

### Individual Alleles

Access through `genotype[index]`:

- `allele.phased`: Boolean indicating if allele is phased
- `allele.allele`: Integer value of allele (0=ref, 1=first alt, etc)

## Sample Objects

Accessed through `variant:sample(sample_name)`:

### Attributes

- All FORMAT fields are available as direct attributes
- `GT`: Special handling for genotype
  - Array of allele values (0-based)
  - `phase`: Array of boolean values indicating phasing
- Common fields (if present in VCF):
  - `DP`: Depth
  - `GQ`: Genotype Quality
  - `AD`: Allelic Depths
  - `PL`: Phred-scaled Likelihoods

## Header Object

Available in prelude scripts for header manipulation:

### Attributes

- `header.samples`: Array of sample names (get/set)

### Methods

- `header:info_get(field_name)`: Get INFO field definition
- `header:format_get(field_name)`: Get FORMAT field definition
- `header:add_info({ID, Number, Type, Description})`: Add new INFO field
- `header:add_format({ID, Number, Type, Description})`: Add new FORMAT field
- `header:add_filter({ID, Description})`: Add new FILTER

### Field Definition Parameters

When adding new fields:

- `ID`: Field identifier
- `Number`: Number of values ("1", "A", "R", "G", etc)
- `Type`: Data type ("Integer", "Float", "String", "Flag")
- `Description`: Field description

## Examples

```lua
-- Filter heterozygous variants
return variant.genotypes[1].alts == 1

-- Get sample depth
local dp = variant:sample("NA12878").DP

-- Add new INFO field in prelude
header:add_info({
    ID="AF_MAX",
    Number="1",
    Type="Float",
    Description="Maximum allele frequency"
})

-- Set complex filter
return variant.qual > 30 and variant:info("DP") > 10 and variant.genotypes[1].alts > 0
```
