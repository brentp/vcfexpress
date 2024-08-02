# load the csq code (csq.lua), then add the af_copy field and define `desc` from the header (pre.lua)
./target/debug/vcfexpress filter -e "csqs = CSQS.new(variant:info('vep'), desc); return csqs:any(function(c) return c.IMPACT == 'HIGH' end)" \
        -o var.bcf gnomad.genomes.v4.0.sites.chrY.vcf.bgz \
        -p scripts/csq.lua \
        -p scripts/pre.lua \
        -s 'af_copy=return variant:info("AF", 0)' 

time ./target/release/vcfexpress filter -e "csqs = CSQS.new(variant:info('vep'), desc); return csqs:any(function(c) return c.IMPACT == 'HIGH' end)" \
    -s "HIGH_IMPACT=csqs = CSQS.new(variant:info('vep'), desc); return csqs:any(function(c) return c.IMPACT == 'HIGH' end)" \
    -p scripts/csq.lua \
    -p scripts/high_impact_header.lua \
    -o hi.bcf gnomad.genomes.v4.0.sites.chrY.vcf.bgz \
