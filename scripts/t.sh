# load the csq code (csq.lua), then add the af_copy field and define `desc` from the header (pre.lua)
./target/debug/vcfexpr filter -e "csq = CSQ.new(variant:info('vep', 0), desc); return csq.IMPACT == 'HIGH'" \
        -o var.bcf gnomad.genomes.v4.0.sites.chrY.vcf.bgz \
        -p scripts/csq.lua \
        -p scripts/pre.lua \
        -s 'af_copy=return variant:info("AF", 0)' 
