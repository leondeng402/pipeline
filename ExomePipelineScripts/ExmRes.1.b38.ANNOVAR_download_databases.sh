#!/bin/bash


dbToDownload=$1
#echo $dbToDownload
AnnovarDirectory=/home/local/users/ld402/resources/annovar/

# 20180801 update
latest1KG=1000g2015aug
latestdbSNP=avsnp150
latestLJB=ljb26_all
latestESPall=esp6500siv2_all
latestESPaa=esp6500siv2_aa
latestESPea=esp6500siv2_ea
latestCOMIC=cosmic70

cd $AnnovarDirectory

case $dbToDownload in
    1)
        nam="Download: refseq hg38 gene reference"
        cmd="annotate_variation.pl -downdb -buildver hg38 -webfrom annovar refGene humandb/"
        ;;
    2)
        nam="Download: 1000 genomes reference"
        cmd="annotate_variation.pl -downdb -buildver hg38 -webfrom annovar $latest1KG humandb/"
        ;;
    3)
        nam="Download: dbSNP reference"
        cmd="annotate_variation.pl -downdb -buildver hg38 -webfrom annovar $latestdbSNP humandb/"
        ;;
    4)
        nam="Download: LJB"
        cmd="annotate_variation.pl -downdb -buildver hg38 -webfrom annovar $latestLJB humandb/"
        ;;
    5)
        nam="Download: ESP alternative allele frequency - all"
        cmd="annotate_variation.pl -downdb -buildver hg38 -webfrom annovar $latestESPall humandb/"
        ;;
    6)
        nam="Download: ESP alternative allele frequency - African Americans"
        cmd="annotate_variation.pl -downdb -buildver hg38 -webfrom annovar $latestESPaa humandb/"
        ;;
    7)
        nam="Download: SP alternative allele frequency - European Americans"
        cmd="annotate_variation.pl -downdb -buildver hg38 -webfrom annovar $latestESPea humandb/"
        ;;
    8)
        nam="Download: superdups "
        cmd="annotate_variation.pl -downdb -buildver hg38 genomicSuperDups humandb/"
        ;;
    9)
        nam="Download: Cadd top 10% "
        cmd="annotate_variation.pl -downdb -buildver hg38 -webfrom annovar caddgt10 humandb/"
        ;;
    10)
        nam="Download: Cadd indel "
        cmd="annotate_variation.pl -downdb -buildver hg38 -webfrom annovar caddindel humandb/"
        ;;
    11)
        nam="Download: exac03"
        cmd="annotate_variation.pl -downdb -buildver hg38 -webfrom annovar  exac03 humandb/"
        ;;
    12)
        nam="Download: Cadd full "
        cmd="annotate_variation.pl -downdb -buildver hg38 -webfrom annovar cadd humandb/"
        ;;
    13)
        nam="Download: clinvar "
        cmd="annotate_variation.pl -downdb -buildver hg38 -webfrom annovar clinvar_20180603 humandb/"
        ;;
    14)
        nam="Download: COSMIC "
        cmd="annotate_variation.pl -downdb -buildver hg38 -webfrom annovar $latestCOMIC humandb/"
        ;;
esac

echo "Start $nam - `date`"
echo $cmd
eval $cmd
if [[ $? == 0 ]]; then
    echo "Finish $nam - `date`"
else
    echo "Error during "${nam/load:/loading}
fi
