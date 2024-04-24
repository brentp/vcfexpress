# vcfexpr

[![Rust](https://github.com/brentp/vcfexpr/actions/workflows/rust.yml/badge.svg)](https://github.com/brentp/vcfexpr/actions/workflows/rust.yml)

This is an experiment on how to implement user-expressions
that can filter (and modify) a VCF and specify an output template.
It uses lua as the expression language. It is fast.

For the optional output template, it uses [luau string templates](https://luau-lang.org/syntax#string-interpolation)
where luau is lua with some extensions and very good speed.

# Examples


extract a single variant and output a bed of the variant:
```
vcfexpr filter -e "return variant.id == 'rs2124717267'" \
    --template '{variant.chrom}\t{variant.start}\t{variant.stop}' -o var.bed $vcf
```
---
filter based on INFO and write bcf:
```
vcfexpr filter -e "return variant:info('AN') > 3000" \
   -o high_an.bcf $input_vcf
```

---
check the sample fields to get variants where `all` samples have high DP.
`all` is defined by `vcfexpr` (`any`, `filter` are also available).
Users can load their own functions with `-l $lua_file`.
```
vcfexpr filter \
   -e 'return all(function (dp)  return dp > 10 end, variant:format("DP"))' \
   -o all-high-dp.bcf $input_vcf
```
---

get all of the FORMAT fields for a single sample into a lua table.
find variant that are high-quality hom-alts.

```
vcfexpr filter \
   -e 's=variant:sample("NA12878"); return s.DP > 10 and s.GQ > 20 and s.GT[1] == 1 and s.GT[2] == 1' \
   -o output.bcf \
   input.vcf
``` 

---

add a new info field (`af_copy`) and set it.
```
$ cat pre.lua
header:add_info({ID="af_copy", Number=1, Description="adding a single field", Type="Float"})
```
then run with:
```
vcfexpr filter -p pre.lua -e 'return variant:format("AD")[1][2] > 0' \
   -s 'af_copy=return variant:info("AF", 0)' \
   input.vcf > output.vcf
```



# Attributes / Functions

```lua
variant.chrom -> string
variant.REF (get/set) -> string
variant.ALT (get/set) -> vec<string>
variant.id (get/set) -> string
variant.start -> integer
variant.stop -> integer
variant.pos (get/set) -> integer -- 0-based
variant.qual (get/set) -> number
variant.filters (get/set) -> vec<string>
variant.FILTER (get/set) -> string (only first one reported)
variant.genotypes -> vec<Genotype>
variant:format("field_name") -> vec<string|number>
-- optional 0-based 2nd arg to info() gets just the desired index.
variant:info("field_name") -> number|string|bool|vec<number|string|bool> 
-- useful to pprint(variant:sample("mysample")) to see available fields.
variant:sample("sample_name") -> table<string=any> 

genotypes = variant.genotypes
genotype = genotypes[i] -- get single genotype for 1 sample
tostring(genotype) -- e.g. "0/1"
genotype.alts -- integer for number of non-zero, non-unknown alleles

allele = genotype[1] 
allele.phased -> bool
allele.allele -> integer e.g. 0 for "0" allele

header.samples (set/get) -> vec<string> -- TODO: allow setting samples before iteration.

-- these header:add_* are available only in the prelude. currently only Number=1 is supported.
header:add_info({Type="Integer", Number=1, Description="asdf", ID="new field"})
header:add_format({Type="Integer", Number=1, Description="xyz", ID="new format field"})


sample = variant:sample("NA12878")
sample.DP -- any fields in the row are available. special case for GT. use pprint to see structure:
pprint(sample)
--[[
{  .GQ = 63,
  .DP = 23,
  .GT = { -- GT gives index into alt alles (or -1 for .)
    [1] = 0,
    [2] = 1},
  .AD = {
    [1] = 23,
    [2] = 0},
  .PL = {
    [1] = 0,
    [2] = 63,
    [3] = 945},
  -- this is the genotype phase.  so with GT, this is 0|1
  .phase = { 
    [1] = false, 
    [2] = true}}
--]]
```


  

# Usage

```
Filter a VCF/BCF and optionally print by template expression. If no template is given the output will be VCF/BCF

Usage: vcfexpr filter [OPTIONS] <PATH>

Arguments:
  <PATH>  Path to input VCF or BCF file

Options:
  -e, --expression <EXPRESSION>
          boolean Lua expression(s) to filter the VCF or BCF file
  -s, --set-expression <SET_EXPRESSION>
          expression(s) to set existing INFO fields (new ones can be added in prelude)
          e.g. --set-expression "AFmax=math.max(variant:info('AF'), variant:info('AFx'))"
  -t, --template <TEMPLATE>
          template expression in luau: https://luau-lang.org/syntax#string-interpolation. e.g. '{variant.chrom}:{variant.pos}'
  -l, --lua <LUA>
          File(s) containing lua code to load. May contain functions that will be called by the expressions
  -p, --lua-prelude <LUA_PRELUDE>
          File containing lua code to run once before any variants are processed
  -o, --output <OUTPUT>
          Optional output file. Default is stdout
  -h, --help
          Print help
```


# TODO

+ Currently --set-expressions can only be used when output is VCF. Update to support template output as well. So we need the header to translate.
+ suuport --set-expressions for FORMAT fields (the infrastructure for this is there, just have to expose it)
+ add a functional lib such as [Moses](https://github.com/Yonaba/Moses) or [Lume](https://github.com/rxi/lume) which have `map`/`filter` and other functions.
  (The user can add these on their own with `--lua`).
+ write a class to simplify accessing CSQ fields.
