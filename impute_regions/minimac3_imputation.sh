#!/bin/bash


counter=0
for i in ${chr[@]}
do
  Minimac3-omp --refHaps ALL.chr${chr[$counter]}.phase3_v5.shapeit2_mvncall_integrated.noSingleton.genotypes.vcf.gz \
          --haps ${prefix[$counter]}_${data}_filtered_phased.vcf \
          --rsid \
          --format GT,DS,GP \
          --chr ${chr[$counter]} \
          --start ${start[$counter]} \
          --end ${end[$counter]} \
          --prefix ${prefix[$counter]}_${data}_imputed \
          --cpus $cpus
        
  counter=$counter+1
done  


. ./postimputation.sh