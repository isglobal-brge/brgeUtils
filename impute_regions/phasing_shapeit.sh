#!/bin/bash

counter=0

for i in ${chr[@]}
do
  shapeit -B $data_chr${chr[$counter]}_filtered \
        -M genetic_map_chr${chr[$counter]}_combined_b37.txt \
        -O $data_chr${chr[$counter]}_filtered_phased \
        --thread $cpus
        
  shapeit -convert \
        --input-haps $data_chr${chr[$counter]}_filtered_phased \
        --output-vcf $data_chr${chr[$counter]}_filtered_phased.vcf
        
  counter=$counter+1
done 


################ I SHOULD HAVE TWO OPTIONS: >200 AND < 200 SAMPLES

# TO PHASE A SMALL NUMBER (<200) OF GWAS SAMPLES
## The following step splits out variants mis-aligned between the reference and gwas panel
#shapeit -check \
#        -V Gwas.chr20.Unphased.vcf\
#        -M genetic_map_chr20.txt \
#        --input-ref reference.haplotypes.gz reference.legend.gz reference.sample \
#        --output-log gwas.alignments

## The following step phases gwas panel using the reference panel while excluding the markers found in the step above.
#shapeit -B gwas \
#        -V Gwas.chr20.Unphased.vcf \
#        --input-ref reference.haplotypes.gz reference.legend.gz reference.sample \
#        --exclude-snp gwas.alignments.strand.exclude \
#        -O Gwas.Chr20.Phased.Output