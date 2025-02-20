# v0.3.5

+ add logging of number of variants/second and number of variants written.
+ add example of groups
+ update to mlua 0.10.3

# v0.3.4

+ Add floats to sample extracted from `variant:sample("MY-SAMPLE")`
+ Expose read-only header during variant evaluation
+ add `samples = variant:samples()` to get all samples and then use as, e.g.

  ```lua

samples = variant:samples()
samples.NA12878.alts == 1 and samples.NA12879.alts == 0 and samples.NA12878.DP > 20
-- to only include some fields. here fields other than `DP` and `GT` are ignored
samples = variant:samples({DP=true})

```
