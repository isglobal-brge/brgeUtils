#!/bin/bash


counter=0

for i in ${chr[@]}
do
  Minimac3-omp --refHaps ALL.chr${chr[$counter]}.phase3_v5.shapeit2_mvncall_integrated.noSingleton.genotypes.vcf.gz \
          --haps ${chr[$counter]}_${data}_filtered_phased.vcf \
          --rsid \
          --format GT,DS,GP \
          --chr ${chr[$counter]} \
          --prefix ${chr[$counter]}_${data}_imputed \
          --cpus $cpus
        
  counter=$counter+1
done  

echo -e "\n\nImputation finished\n\n"