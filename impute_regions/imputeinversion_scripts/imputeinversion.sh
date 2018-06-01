#!/bin/bash

usage() { echo -e "\nUsage: $0 \t[-f,--dataformat <'PLINK' or 'VCF'>]\n\t\t\t\t[-d,--data <Files name without extension (format must be PLINK or VCF)>]\n\t\t\t\t[-i,--inversions <File containing inversion ranges | Inversion range | Inversion name | 'All'>]\n\t\t\t\t[-t,--threads <Integer. Number of threads to be used while phasing with SHAPEIT and while imputing with Minimac3>]\n\t\t\t\t[-k,--keep-files <Intermediate files created during the process will be removed by default. In order to keep them, indicate it writing 'Yes', 'YES' or 'yes'>]\n\nValid inversion names: ${SORTEDKEYS[@]}\n" 1>&2; exit 1; }

declare -A RANGES=( ["inv8p23.1"]="8:8055789-11980649" ["inv17q21.31"]="17:43661775-44372665" ["inv7p11.2"]="7:54290974-54386821" ["invXq13.2"]="X:72215927-72306774" 
                    ["inv12_004"]="12:47290470-47309756" ["inv7_011"]="7:70426185-70438879" ["inv7_003"]="7:31586765-31592019" ["inv11_001"]="11:41162296-41167044" 
                    ["inv2_013"]="2:139004949-139009203" ["inv6_006"]="6:130848198-130852318" ["inv3_003"]="3:162545362-162547641" ["inv7_014"]="7:151010030-151012107" 
                    ["inv11_004"]="11:66018563-66019946" ["inv1_008"]="1:197756784-197757982" ["inv16_017"]="16:85188639-85189823" ["inv21_005"]="21:28020653-28021711" 
                    ["inv12_006"]="12:71532784-71533816" ["inv6_002"]="6:31009222-31010095" ["inv14_005"]="14:65842304-65843165" ["inv1_004"]="1:92131841-92132615" 
                    ["inv2_002"]="2:33764554-33765272" )
                    
declare -A REVRANGES=( ["8:8055789-11980649"]="inv8p23.1" ["17:43661775-44372665"]="inv17q21.31" ["7:54290974-54386821"]="inv7p11.2" ["X:72215927-72306774"]="invXq13.2" 
                       ["12:47290470-47309756"]="inv12_004" ["7:70426185-70438879"]="inv7_011" ["7:31586765-31592019"]="inv7_003" ["11:41162296-41167044"]="inv11_001" 
                       ["2:139004949-139009203"]="inv2_013" ["6:130848198-130852318"]="inv6_006" ["3:162545362-162547641"]="inv3_003" ["7:151010030-151012107"]="inv7_014" 
                       ["11:66018563-66019946"]="inv11_004" ["1:197756784-197757982"]="inv1_008" ["16:85188639-85189823"]="inv16_017" ["21:28020653-28021711"]="inv21_005" 
                       ["12:71532784-71533816"]="inv12_006" ["6:31009222-31010095"]="inv6_002" ["14:65842304-65843165"]="inv14_005" ["1:92131841-92132615"]="inv1_004" 
                       ["2:33764554-33765272"]="inv2_002" )
                    
KEYS=(${!RANGES[@]})
IFS=$'\n' SORTEDKEYS=($(sort -t v -k 2 -g<<<"${KEYS[*]}"))
unset IFS


while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -f|--dataformat) ################ Should do something if it is not valid
    format="$2"
    shift # past argument
    shift # past value
    ;;
    -d|--data)
    data="$2"
    shift # past argument
    shift # past value
    ;;
    -i|--inversions)
    inversion="$2"
    shift # past argument
    shift # past value
    ;;
    -t|--threads)
    cpus="$2"
    shift # past argument
    shift # past value
    ;;
    -k|--keep-files)
    keep_files="$2"
    shift # past argument
    shift # past value
    ;;
    -w|--whole-chr) ############################## modify the script to be able to impute the whole chromosome or just the inversion region? (PUT ALL THE OPTIONS IN THE USAGE MESSAGE)
    impute_whole_chr="$2"
    shift # past argument
    shift # past value
    ;;
    -h|--help|*)
    usage
    shift # past argument
    ;;
esac
done

## Exit if no data provided or data format is not valid
if [[ -z "$data" ]]
then
  echo -e "\nNo input data to be imputed"
  usage
elif [[ -z "$format" ]]
then
  echo -e "\nData format not specified"
  usage
elif [[ $format != plink ]] && [[ $format != PLINK ]] && [[ $format != vcf ]] && [[ $format != VCF ]]
then
  echo -e "\nIndicate valid data format"
  usage
fi

## Determining inversions to impute
if [[ -z "$inversion" ]] # If $inversion is empty (no input indicating inversion to impute)
then
  echo -e "\nNo input indicating inversion(s) region(s)"
  usage
elif [[ $inversion == *.txt ]] # If input is a text file containing the inversions
then 
  chr=( $(awk -F '[ \t:-]+' '{print $1}' $inversion) )
  start=( $(awk -F '[ \t:-]+' '{print $2}' $inversion) )
  end=( $(awk -F '[ \t:-]+' '{print $3}' $inversion) )
  counter=0
  getnames=()
  for i in ${chr[@]}
  do
    getnames+=( "${chr[$counter]}:${start[$counter]}-${end[$counter]}" )
    counter=$counter+1
  done
elif [[ " ${RANGES[*]} " == *" ${RANGES[$inversion]} "* ]] # If input is the name of the inversion
then
  chr=${RANGES[$inversion]%:*}
  start=${RANGES[$inversion]##*:}
  start=${start%-*}
  end=${RANGES[$inversion]##*-} 
  prefix=$inversion
elif [[ " ${RANGES[*]} " == *" $inversion "* ]] # If input is the genomic range of the inversion
then
  chr=${inversion%:*}
  start=${inversion##*:}
  start=${start%-*}
  end=${inversion##*-}
  for i in "${!RANGES[@]}"
  do
    if [[ " ${RANGES[$i]} " == " $inversion " ]]
    then
      prefix=$i
    fi
  done
elif [[ $inversion == all ]] || [[ $inversion == All ]] || [[ $inversion == ALL ]] # If input is ALL|All|all
then
  chr=()
  start=()
  end=()
  for i in ${RANGES[@]}
  do
    chr+=( ${i%:*} )
    middle=${i##*:}
    start+=( ${middle%-*} )
    end+=( ${i##*-} )
  done
  counter=0
  getnames=()
  for i in ${chr[@]}
  do
    getnames+=( "${chr[$counter]}:${start[$counter]}-${end[$counter]}" )
    counter=$counter+1
  done
else # If the input does not correspond with any of the previous options, it is not valid 
  echo -e "\nInput indicating inversion(s) region(s) is NOT valid"
  usage
fi

if [[ -z "$prefix" ]]
then
  prefix=()
  for i in "${getnames[@]}"
  do
    if [[ " ${REVRANGES[*]} " == *" ${REVRANGES[$i]} "* ]]  
    then
      prefix+=( ${REVRANGES[$i]} )
    fi
  done
fi

unique_chr=( $(echo "${chr[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' ') )


## Number of threads (By default: $cpus = 1)
if [[ -z "$cpus" ]]
then
  echo -e "\nNumber of threads was not indicated. Value set to 1"
  cpus=1
elif ! [[ "$cpus" =~ ^[0-9]+$ ]]
then
  echo -e "\nInput for number of threads not valid. Value set to 1"
  cpus=1
fi



. ./prephasing.sh





























