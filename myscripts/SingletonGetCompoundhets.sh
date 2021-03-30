while getopts v:t:n:a:p: opt; do
    case "$opt" in
        v) VcfFil="$OPTARG";;
        t) SingFil="$OPTARG";;
        p) Proband="$OPTARG";;
        n) DirPre="$OPTARG";;
        a) ArrNum="$OPTARG";;
        #H) echo "$usage"; exit;;
    esac
done
R --vanilla <<RSCRIPT
options(stringsAsFactors=F)

het <- read.delim("$Proband.tempcHet.tsv")

gens <- unique(het[duplicated(het[,"Gene"]),"Gene"])

comphet <- het[het[,"Gene"]%in%gens,]

write.table(comphet, "$Proband.compound_heterozygous.tsv", col.names=T, row.names=F, quote=F, sep="\t")
RSCRIPT
LEN=`cat $Proband.compound_heterozygous.tsv | wc -l`
if [[ $LEN -gt 1 ]]; then
    nohup $FiltScrDir/xAnnotateVariantTSV.sh -i $Proband.compound_heterozygous.tsv &
fi
mv $Proband.tempcHet.log $Proband.compound_heterozygous.log
rm $Proband.tempcHet*

