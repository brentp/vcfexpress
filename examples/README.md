# Examples

This contains a set of examples of how to use vcfexpress.

We welcome additional contributions!

<details>
<summary>Test for Hardy-Weinberg Equilibrium</summary>

Here we filter out variants that are out of HWE (have a low p-value)

```
vcfexpress filter -p examples/hwe.lua -e "hwe(variant) > 0.05" $vcf -o variants-in-hwe.bcf
```

</details>

<details>
<summary>Filter variants on calculated per-sample posterior genotype probability</summary>

We can use genotype likelihoods to get a posterior probability for each genotype -- and therefore for the
called genotype. We then filter to variants with a high (0.95+) probability of the called genotype.
Here we require only that `any` samples meet that threshold.

The relevant lua code in the functions loaded is:

```lua
    --- ... see full code in examples/gl_prob.lua
    -- Compute sum of all likelihoods
    local sum_L = 0
    for i = 1, #GL do
        sum_L = sum_L + 10 ^ GL[i]
    end

    -- Compute probability for the requested index
    return (10 ^ GL[alts + 1]) / sum_L

```

Which gives scaled probabilities that sum to 1 across possible genotypes and returns the probability for the called genotype.

```
vcfexpress filter \
    -p examples/gl_prob.lua \
    -e 'return any(function(sm) return GT_prob(sm) > 0.95 end, variant:samples({GL=true}))' \
    $vcf
```

Note that in cases like this where we access all samples for a given variant, `variant:samples()` will be much faster
than calling `variant:sample()` for each sample iteratively. Using `variant:format('GL')` can be even faster, but returns
a table (of tables) so is less convenient.

</details>

<details>
<summary>Set the ID field of a variant</summary>

we can set the ID field of the variant. here we use the following lua code in `examples/set-id-to-chrom-start-ref-alt.lua` to do so:

```lua
function set_id(variant)
    local alt = variant.ALT[1]
    variant.id = string.format("%s-%d-%s-%s", variant.chrom, variant.pos, variant.REF, alt)
    return true
end
```

Then we call as:

```
vcfexpress filter \
    -p examples/set-id-to-chrom-start-ref-alt.lua \
    -e "return set_id(variant)" \
    examples/trio.vcf.gz -o set.vcf.gz
```

and the output looks like:

```
$ zcat set.vcf.gz | grep -v ^## | cut -f 1-5 | head
#CHROM  POS     ID      REF     ALT
1       876499  1-876498-A-G    A       G
1       887560  1-887559-A-C    A       C
1       887801  1-887800-A-G    A       G
1       888639  1-888638-T-C    T       C
1       888659  1-888658-T-C    T       C
1       897325  1-897324-G-C    G       C
1       906272  1-906271-A-C    A       C
1       908823  1-908822-G-A    G       A
1       909238  1-909237-G-C    G       C
```

</details>

<details>
<summary>Filter output based on trio genotypes</summary>

Here, we extract sites where all samples are heterozygotes (with 1 alternate allele).

```
vcfexpress filter \
    -e "local gts = variant.genotypes; return gts[1].alts == 1 and gts[2].alts == 1 and gts[3].alts == 1" \
    examples/trio.vcf.gz \
   | grep -v ^## | cut -f 9- | head
GT:DP:RO:AO     0/1:41:15:26    0/1:42:20:22    0/1:53:25:27
GT:DP:RO:AO     0/1:53:28:25    0/1:65:34:31    0/1:57:24:33
GT:DP:RO:AO     0/1:43:22:19    0/1:50:19:31    0/1:49:23:25
GT:DP:RO:AO     0/1:43:18:25    0/1:45:21:24    0/1:57:50:7
GT:DP:RO:AO     0/1:90:63:27    0/1:97:38:59    0/1:100:77:23
GT:DP:RO:AO     0/1:54:43:11    0/1:63:52:11    0/1:71:57:14
GT:DP:RO:AO     0/1:62:28:33    0/1:63:31:32    0/1:56:31:25
GT:DP:RO:AO     0/1:50:22:28    0/1:40:19:21    0/1:56:28:28
GT:DP:RO:AO     0/1:44:21:23    0/1:55:24:31    0/1:38:13:25
GT:DP:RO:AO     0/1:54:25:29    0/1:64:29:35    0/1:39:19:19
```

</details>

<details>
<summary>update FILTER based on INFO</summary>

Here we update the FILTER field based on the variant QUAL field.

First, we add the filter to the header using code run in the prelude:

```
echo 'header:add_filter({ID="LowQual", Description="Qual less than 1000"})' > examples/add_filter_to_header.lua
```

Then we use that and the expression to update the FILTER where appropriate:

```
vcfexpress filter -p examples/add_filter_to_header.lua -e "if variant.qual < 1000 then variant.FILTER = 'LowQual' end; return true" examples/trio.vcf.gz | grep -v ^## | head | cut -f 1-7
#CHROM  POS     ID      REF     ALT     QUAL    FILTER
1       876499  .       A       G       18274.6 .
1       887560  .       A       C       28966.1 .
1       887801  .       A       G       25116.2 .
1       888639  .       T       C       24827.5 .
1       888659  .       T       C       23801   .
1       897325  .       G       C       23174.8 .
1       906272  .       A       C       5300.31 .
1       908823  .       G       A       875.757 LowQual
1       909238  .       G       C       6998.41 .
```

The output header also contains:

```
##FILTER=<ID=LowQual,Description="Qual less than 1000">
```

</details>
