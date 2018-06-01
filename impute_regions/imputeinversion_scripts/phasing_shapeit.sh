#!/bin/bash



#numberind=$(wc -l $data.fam | cut -d' ' -f1)
#if [[ $numberind > 200 ]]
#then
#  echo "The file has more than 200 individuals"
#  echo "$numberind"
#else
#  echo "The file has NO more than 200 individuals"
#  echo "$numberind"
#fi

mkdir phased_files

for i in ${unique_chr[@]}
do
  mkdir phased_files/${i}
  if [[ $i == X ]]
  then
    shapeit -B prephasing_files/${i}/${i}_${data}_filtered \
          -M /home/isglobal.lan/itolosana/homews/reference_panels/genetic_map_chrX_nonPAR_combined_b37.txt \
          -O phased_files/${i}/${i}_${data}_filtered_phased \
          --chrX \
          --thread $cpus
          
    shapeit -convert \
          --input-haps phased_files/${i}/${i}_${data}_filtered_phased \
          --output-vcf phased_files/${i}/${i}_${data}_filtered_phased.vcf
  else
    shapeit -B prephasing_files/${i}/${i}_${data}_filtered \
          -M /home/isglobal.lan/itolosana/homews/reference_panels/genetic_map_chr${i}_combined_b37.txt \
          -O phased_files/${i}/${i}_${data}_filtered_phased \
          --thread $cpus
        
    shapeit -convert \
          --input-haps phased_files/${i}/${i}_${data}_filtered_phased \
          --output-vcf phased_files/${i}/${i}_${data}_filtered_phased.vcf
  fi
done 

mv shapeit_* phased_files/

if [[ -z "$keep_files" ]] || [[ $keep_files != Yes ]] && [[ $keep_files != YES ]] && [[ $keep_files != yes ]] # By default, the intermediate files will be deleted
then
  rm -rf prephasing_files
fi

. ./minimac3_imputation.sh
################################################################## ELIMINAR LOS ARCHIVOS QUE GENERA SHAPEIT, LOS .LOG, LOS .FRQ Y LOS .HH

#shapeit --input-vcf gwas.vcf \
#        -M genetic_map.txt \
#        -O gwas.phased

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


