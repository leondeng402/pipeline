#!/bin/bash

# set default arguments
Usage="illimina_filter.sh -i <IlluminaVCFFile> -o <OutputFile>

       -i: illumina vcf file
       -o: output file
"

while getopts i:o: opt; do
    case "$opt" in
        i) in_file="$OPTARG";;
        o) out_file="$OPTARG";;
    esac
done

echo "input = $in_file" 
echo "output = $out_file"

line_num=0
while read -r line
do
    ((line_num++))

    #Write the vcf header
    #If line starts with #, write to the outfile
    if [[ $line == \#* ]]; then
        if [[ $line_num -eq 1 ]]; then
            echo "$line" > $out_file
        else
            echo "$line" >> $out_file
        fi
        continue
    fi
    
    #Check Field7_FILTER, Field10-**_samples-genotypes    
    field_num=0
    line_check=0
    #Check the FILTER==PASS, Genotype!=./.
    for i in $(echo $line | tr "\t" "\n")
        do
            ((field_num++))
            # Check Filter whether it is a PASS
            if [[ $field_num -eq 7 ]]; then
                if [[ $i == "PASS" ]]; then
                    line_check=1
                else
                    line_check=0
                fi      
            fi

            # Check samples genotypes
            if [[ $field_num -ge 10 ]] && [[ $line_check -eq 1 ]]; then
                element_num=0
                for j in $(echo $i | tr ":" "\n")
                    do
                        ((element_num++))
                        if [[ $element_num -eq 1 ]]; then
                           # check the genotype format
                           if [[ $j == *[/]* ]]; then
                              line_check=1
                           else
                              line_check=0
                              echo "$j"
                           fi

                           if [[ $j == './.' ]]; then
                               line_check=0
                           fi
                           
                        fi
                    done
            fi
           
           # Check the genotype quality
          #if [[ $field_num -eq 9 ]] && [[ $line_check -eq 1 ]]; then 
          #    echo "check the genotype qualities, to be implemented"
          #fi
        done

    if [[ $line_check -eq 1 ]]; then
        echo "$line" >> $out_file
    fi
    #echo "Filter the quality" 
    #echo "Filter the genotype"
   
done < "$in_file"
