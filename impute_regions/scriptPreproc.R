#'#################################################################################
# Example of preprocsesing data for imputation ####
#'#################################################################################
plink --bfile bladder_all --recode --out bladder

sed 's/^23/X/' bladder.map  > bladder18.map 
cp bladder.ped bladder18.ped

# Convert hg18 to hg19
python /home/isglobal.lan/cruiz/liftOverPlink/liftOverPlink.py --map bladder18.map --out lifted --chain ~/liftOverPlink/hg18ToHg19.over.chain.gz
python /home/isglobal.lan/cruiz/liftOverPlink/rmBadLifts.py --map lifted.map --out good_lifted.map --log bad_lifted.dat
cut -f 2 bad_lifted.dat > to_exclude.dat
cut -f 4 lifted.bed.unlifted | sed "/^#/d" >> to_exclude.dat 
plink --file bladder18 --recode ped --out lifted --exclude to_exclude.dat 
plink --ped lifted.ped --map good_lifted.map --make-bed --out final

# Remove samples with low call rate
plink --bfile final --mind 0.1 --make-bed --out final_filtered


# Create files for imputation
plink --bfile final --freq --out final
perl ~/HRC-1000G-check-bim-v4.2.pl -b final.bim -f final.frq -r ../../HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h

## Execute first lines of Run_plink.sh
## Create files for imputation
array=(1 2 3 6 7 8 16 17 21)
for i in "${array[@]}"
do
plink --bfile final-updated --chr $i --chr-output M --set-hh-missing --recode vcf-iid --out Bladder-chr$i
vcf-sort Bladder-chr$i.vcf | bgzip -c > Bladder-chr$i.recode.vcf.gz
done


## Create files for peddy
plink --bfile final --recode vcf-iid bgz --out merged
tabix -p vcf merged.vcf.gz

## Run peddy
python -m peddy --prefix cancer merged.vcf.gz peddy.fam
