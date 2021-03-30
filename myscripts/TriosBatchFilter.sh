#!/bin/bash
# set default arguments
Usage="TrisoBatchFilter.sh -i <TrioListFile> -v <VariantsFile> 

	-i (required): the file containing the list of trios
	-v (required): the compressed variants vcf file in the format of .vcf.gz
        -D (flag): switch on the denovo quality
	-H (flag): echo the usage message and exit

"

DeNovoQual=false
while getopts i:v:n:HD opt; do 
    case "$opt" in 
        i) TrioListFile="$OPTARG";;
        v) VariantsFile="$OPTARG";;
        D) DeNovoQual=true;;
	H) echo "$Usage"; exit;;
    esac
done
#echo "$TrioListFile"
#echo "$VariantsFile"

TriosNum=$(cat $TrioListFile | wc -l);
#echo "$TriosNum"
echo "DeNovoQual=$DeNovoQual"

# set the filtering scripts directory
FilterScriptsDir="/home/local/users/ld402/scripts/Filtering_scripts"

# check the input
if [[ ! -e "$TrioListFile" ]]; 
then
    echo "Error: Missing required arguments -i"; echo "$Usage"; exit;
fi

if  [[ ! -e "$VariantsFile" ]]; 
then    
    echo "Error: Missing required arguments -v"; echo "$Usage"; exit;
fi
   
if  [[ $TriosNum -lt 1 ]]; 
then    
    echo "Error: Trios number is less than 1";  exit;
fi

# run trios filter for each trio
for (( i=1; i <= $TriosNum; i++))
do
    #echo "i=$i"
    if [ "$DeNovoQual" == "true" ]; then
        cmd="$FilterScriptsDir/xTrioFilters.sh -v $VariantsFile -t $TrioListFile -a $i -D"
    else
        cmd="$FilterScriptsDir/xTrioFilters.sh -v $VariantsFile -t $TrioListFile -a $i"
    fi
    echo "command: $cmd"
    eval $cmd
done

