#!/bin/bash

# set default arguments
Usage="DOCBatchCompute.sh -i <FileList> -r <ReferenceFile> -t <TargetFile> -l <LogFile> -BFH

        -i (required): the file containing the list of paths for the bam files
	-r (required): shell command with the path for reference and resource
	-t (required): exome capture kit targets
	-l (optional): log file
	-B (flag): prevent GATK from phoning home
	-F (flag): Fix mis-encoded base quality scores
	-H (flag): echo the usage message and exit
"
BadET="false"
FixMisencoded="false"
while getopts i:r:l:t:BFH opt; do
    case "$opt" in 
        i) InputFile="$OPTARG";;
        r) RefFile="$OPTARG";;
        t) TgtBed="$OPTARG";;
	l) LogFile="$OPTARG";;
        B) BadET="true";;
        F) FixMisencoded="true";;
        H) echo "$Usage"; exit;;
    esac
done
pipe_dir="/home/local/users/ld402/scripts/ExomePipelineScripts" 
# check the required variable
if [[ ! -e "$InputFile" ]] || [[ ! -e "$RefFile" ]] || [[ -e "$TgtBed" ]]; 
then
    echo "Missing required arguments"; echo "$Usage"; exit;
fi

for i in $(cat $InputFile); do
    bam_base=${i##*/}
    echo $bam_base
    StepCmd="$pipe_dir/ExmAln.8a.DepthofCoverage.sh -i $i -r $RefFile -t $TgtBed -BD"
    printf "$StepCmd\n\n"
    
    eval $StepCmd
done



    
