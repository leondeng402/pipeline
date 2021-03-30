#!/bin/bash

#This script initiates the bam-->fastq-->bam-->gVCF realignment pipeline from a file containing a list of bam files.

#    InpFil - (required) - Path to Bam file to be aligned. Alternatively a file with a list of bams can be provided and the script run in an array job. List file name must end ".list"
#    RefFil - (required) - shell file containing variables with locations of reference files and resource directories; see list below for required variables
#    LogFil - (optional) - File for logging progress
#    TgtBed - (optional) - Exome capture kit targets bed file (must end .bed for GATK compatability) ; may be specified using a code corresponding to a variable in the RefFil giving the path to the target file- only required if calling pipeline
#    PipeLine - P -(flag) - will start the GATK realign and recalibrate pipeline using the files generated by this script
#    Flag - F - Fix mis-encoded base quality scores - Rewrite the qulaity scores from Illumina 1.5+ to Illumina 1.8+, will subtract 31 from all quality scores; used to fix encoding in some datasets (especially older Illumina ones) which starts at Q64 (see https://en.wikipedia.org/wiki/FASTQ_format#Encoding)
#    Help - H - (flag) - get usage information


#list of variables required in reference file:
# $EXOMPPLN - directory containing exome analysis pipeline scripts, 

#set default arguments
usage="
ExmAln.0b.StartBamRealingmentFromList.sh -i <InputFile> -r <reference_file> -t <target intervals file> -l <logfile> -PRFH

     -i (required) - List of Bam files, or 3 column table for starting Fastq alignment - see ExmAln.1a.Align_Fastq_to_Bam_with_BWAmem.sh for details
     -r (required) - shell file containing variables with locations of reference files and resource directories
     -l (optional) - Log file
     -t (optional) - Exome capture kit targets or other genomic intervals bed file (must end .bed for GATK compatability); this file is required if calling the pipeline but otherwise can be omitted
     -P (flag) - Initiate exome analysis pipeline after completion of script; otherwise stop at realigned bam
     -F (flag) - Fix mis-encoded base quality scores - Rewrite the qulaity scores from Illumina 1.5+ to Illumina 1.8+
     -H (flag) - echo this message and exit
"

#get arguments
while getopts i:r:l:t:PFH opt; do
    case "$opt" in
        i) InpFil="$OPTARG";;
        r) RefFil="$OPTARG";; 
        l) LogFil="$OPTARG";;
        t) TgtBed="$OPTARG";; 
        P) PipeLine="true";;
        F) FixMisencoded="true";;
        H) echo "$usage"; exit;;
    esac
done

#Call the RefFil to load variables
RefFil=`readlink -f $RefFil`
source $RefFil

#LogFile
LogFil=StartRealignFromList.$$.log

#check input file type - bam or fastq
FirstLine=`head -n 1 $InpFil`
NCOL=`head -n 1 $InpFil | tr '\t' '\n' | wc -l`
echo $NCOL
echo $FirstLine
if [[ $FirstLine == *bam && NCOL -eq 1 ]]; then
    FileType="bam"
elif [[ $NCOL -eq 3 ]]; then
    FirstLine=`head -n 1 $InpFil | cut -f 1`
    echo $FirstLine
    if [[ $FirstLine == *fq* || $FirstLine == *fastq* ]]; then
        FileType="fastq"
    fi
fi
if [[ ! $FileType ]]; then
    echo "Error with input file"
    exit 1
fi

echo "File type is "$FileType
## This file also requires exome.lib.sh - which contains various functions used throughout the Exome analysis scripts; this file should be in the same directory as this script

if [[ $FileType == "bam" ]]; then
    for i in $(cat $InpFil); do
        BamNam=`basename $i`
        StepCmd="nohup $EXOMPPLN/ExmAln.1b.ReAlign_Bam_with_BWAmem.sh -i $i -t $TgtBed -r $RefFil"
        echo "$StepCmd"
        if [[ $PipeLine == "true" ]]; then StepCmd=$StepCmd" -P"; fi
        if [[ $FixMisencoded == "true" ]]; then StepCmd=$StepCmd" -F"; fi
        StepCmd=$StepCmd" > RealignBam.$BamNam.o 2>&1 &"
        echo $StepCmd >> $LogFil
        eval $StepCmd
    done
else
    FILLEN=`cat $InpFil | wc -l`
    for i in $(seq 1 $FILLEN); do
        FilNam=$(head -n $i $InpFil | tail -n 1 | cut -f 2 | sed s/.*SM:// | sed 's/\\t.*//')
        StepCmd="nohup $EXOMPPLN/ExmAln.1a.Align_Fastq_to_Bam_with_BWAmem.sh -i $InpFil -a $i -t $TgtBed -r $RefFil"
        echo "$StepCmd"
        if [[ $PipeLine == "true" ]]; then StepCmd=$StepCmd" -P"; fi
        if [[ $FixMisencoded == "true" ]]; then StepCmd=$StepCmd" -F"; fi
        StepCmd=$StepCmd" > AlignFastq.$FilNam.o 2>&1 &"
        echo $StepCmd >> $LogFil
        eval $StepCmd
    done
fi
