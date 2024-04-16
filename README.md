# vcfexpr

[![Rust](https://github.com/brentp/vcfexpr/actions/workflows/rust.yml/badge.svg)](https://github.com/brentp/vcfexpr/actions/workflows/rust.yml)

This is an experiment on how to implement user-expressions
that can filter (and soon modify) a VCF and specify an output template.
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
use `all` in user-defined `funcs.lua`
check the sample fields to get variants where all samples have high DP
```
vcfexpr filter -l funcs.lua \
   -e 'return all(variant:format("DP"), function (dp)  return dp > 10 end)' \
   -o all-high-dp.bcf $input_vcf
```
```
$ cat funcs.lua
function all(t, f)
    for _, v in pairs(t) do
        if v ~= nil and not f(v) then
            return false
        end
    end
return true
```

# Usage

```
Filter the VCF or BCF file and optionally apply a template. If no template is given the output will be VCF/BCF (TODO)

Usage: vcfexpr filter [OPTIONS] <PATH>

Arguments:
  <PATH>  Path to input VCF or BCF file

Options:
  -e, --expression <EXPRESSION>  boolean Lua expression(s) to filter the VCF or BCF file
  -t, --template <TEMPLATE>      template expression in luau: https://luau-lang.org/syntax#string-interpolation. e.g. '{variant.chrom}:{variant.pos}'
  -l, --lua <LUA>                File(s) containing lua code to run. Can contain functions that will be called by the expressions
  -o, --output <OUTPUT>          optional output file. Default is stdout
  -h, --help                     Print help
```
