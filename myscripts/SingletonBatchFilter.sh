#!/bin/bash

# set default arguments
Usage="SingletonBatchFilter.sh -v <VCFFile> -t <TableFile> -n <Number>

       -v: vcf file
       -t: file containing singletons with fields of familyID, singletonID
       -n: total number of singletons
       -h: echo this help message    
"

while getopts v:t:n:h opt; do
    case "$opt" in
        v) VCFFile="$OPTARG";;
        t) TableFile="$OPTARG";;
        n) NumOfSingles="$OPTARG";;
        h) echo "$Usage"; exit;;
    esac
done
printf "$VCFFile, $TableFile, $NumberOfSingles"

for i in $(seq 1 "$NumOfSingles"); do
    StepCmd="/home/local/users/ld402/scripts/Filtering_scripts/xSingletonFilters.sh -v $VCFFile -t $TableFile -a $i"
    printf "$StepCmd\n\n"
    eval $StepCmd
done
