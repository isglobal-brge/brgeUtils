#!/bin/bash

declare -A RANGES=( ["inv8p23.1"]="8:8055789-11980649" ["inv17q21.31"]="17:43661775-44372665" ["inv7p11.2"]="7:54290974-54386821" ["invXq13.2"]="X:72215927-72306774" 
                    ["inv12_004"]="12:47290470-47309756" ["inv7_011"]="7:70426185-70438879" ["inv7_003"]="7:31586765-31592019" ["inv11_001"]="11:41162296-41167044" 
                    ["inv2_013"]="2:139004949-139009203" ["inv6_006"]="6:130848198-130852318" ["inv3_003"]="3:162545362-162547641" ["inv7_014"]="7:151010030-151012107" 
                    ["inv11_004"]="11:66018563-66019946" ["inv1_008"]="1:197756784-197757982" ["inv16_017"]="16:85188639-85189823" ["inv21_005"]="21:28020653-28021711" 
                    ["inv12_006"]="12:71532784-71533816" ["inv6_002"]="6:31009222-31010095" ["inv14_005"]="14:65842304-65843165" ["inv1_004"]="1:92131841-92132615" 
                    ["inv2_002"]="2:33764554-33765272")
KEYS=(${!RANGES[@]})


while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -f|--dataformat)
    format="$2"
    shift # past argument
    shift # past value
    ;;
    -d|--data) ############################# IF DATA IS EMPTY: ALGÚN ECHO Y USAGEEEEE
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
    -o|--output-prefix)
    prefix="$2"
    shift # past argument
    shift # past value
    ;;
    -h|--help|*)
    usage
    shift # past argument
    ;;
esac
done


if [[ -z "$inversion" ]] # If $inversion is empty (no input indicating inversion to impute)
then
  echo -e "\nNo input indicating inversion(s) region(s)"
  usage
elif [[ $inversion == *.txt ]] # If input is a text file containing the inversions
then 
  chr=( $(awk -F '[ \t:-]+' '{print $1}' $inversion) )
  start=( $(awk -F '[ \t:-]+' '{print $2}' $inversion) )
  end=( $(awk -F '[ \t:-]+' '{print $3}' $inversion) )
elif [[ " ${RANGES[*]} " == *" ${RANGES[$inversion]} "* ]] # If input is the name of the inversion
then
  chr=${RANGES[$inversion]%:*}
  start=${RANGES[$inversion]##*:}
  start=${start%-*}
  end=${RANGES[$inversion]##*-} 
elif [[ " ${RANGES[*]} " == *" $inversion "* ]] # If input is the genomic range of the inversion
then
  chr=${inversion%:*}
  start=${inversion##*:}
  start=${start%-*}
  end=${inversion##*-} 
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
else # If the input does not correspond with any of the previous options, it is not valid 
  echo -e "\nInput indicating inversion(s) region(s) is NOT valid"
  usage
fi


# I guess it is not the most important thing or really relevant for what it will take me, but it is nice and the good thing would be to implement this also (mainly, actually) when there is a list of inversions instead of just one inversion
somet="${chr[@]}:${start[@]}-${end[@]}" # I WILL HAVE TO DO THIS WITH A LOOP, BECAUSE I WILL HAVE TO INDICATE THAT IT IS "${chr[1]}:${start[1]}-${end[1]}" (ONE BY ONE)
echo "${somet[@]}"
if [[ -z "$prefix" ]]
then
  prefix=()
  echo "Prefix is empty"
  for i in "${!RANGES[@]}"
  do
    if [[ " ${RANGES[$i]} " == " $somet[*] " ]]
    then
      prefix+=( $i )
    fi
  done
  #prefix=(${!RANGES[$chr:$start-$end]})
  #echo "${RANGES[2:33764554-33765272]}"
  #prefix="${chr[$counter]}_imputed_final.vcf.gz" # Algo habrá que hacer aquí porque si piden all, serán los mismos cromosomas (ponerlo con el nombre de la inversión? yasyasssss -> cómoooooo?? )
else
  echo "I said my files are $prefix"
fi

echo "Somehow it works, ${prefix[@]}"












