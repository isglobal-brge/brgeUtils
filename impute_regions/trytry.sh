#!/bin/bash


while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -f|--dataformat)
    newformat="$2"
    shift # past argument
    shift # past value
    ;;
    -d|--data)
    newdata="$2"
    shift # past argument
    shift # past value
    ;;
    -i|--inversions)
    inversion="$2"
    shift # past argument
    shift # past value
    ;;
    -t|--thread)
    cpus="$2"
    shift # past argument
    shift # past value
    ;;
    -h|--help|*)
    usage
    shift # past argument
    ;;
esac
done



if [[ $newformat == plink ]] || [[ $newformat == PLINK ]]
then
  echo "I know $data is PLINK from the other file."
  newvariable="hello there"
elif [[ $newformat == vcf ]] || [[ $newformat == VCF ]]
then
  echo "I know $data is VCF from the other file."
fi 

counter=0
for i in ${chr[@]}
do
  echo "${chr[$counter]}"  
  counter=$counter+1
done  

numberind=$(wc -l $data.fam | cut -d' ' -f1)

if [[ $numberind > 200  ]]
then
  echo "The file has more than 200 individuals"
  echo "$numberind"
else
  echo "The file has NO more than 200 individuals"
  echo "$numberind"
fi



















