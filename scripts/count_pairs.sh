set -e 
set -u

if [ $# -ne 3 ]
then
    echo "$0 <BAM> <GTF> <OUTDIR>"
    exit 1
fi

curfile=$1
annotation=$2
output=$3
prefix=$(basename $curfile | sed "s/.bam$//")
if [ ! -e $prefix ]
then
    mkdir ${prefix}
    StripeThis ${prefix}
fi

# Stripping out supplementary alignments and putting each read in its own file.
b1=${prefix}/1.bam
b2=${prefix}/2.bam
samtools view -h -F 2048 -f 64 $curfile -b -o $b1
samtools view -h -F 2048 -f 128 $curfile -b -o $b2

# Running through featureCounts.
for x in 1 2
do 
    echo "
library(Rsubread)
annotation <- normalizePath('${annotation}')
setwd('${prefix}')
stats <- featureCounts('${x}.bam', annot.ext=annotation,
    isGTFAnnotationFile=TRUE, minMQS=255, reportReads='CORE', ignoreDup=TRUE)
" | R --no-save --vanilla
done

# Checking that the rownames are the same.
tmp1=${prefix}/1.txt
tmp2=${prefix}/2.txt
cut -f 1 ${b1}.featureCounts > $tmp1
cut -f 1 ${b2}.featureCounts > $tmp2
if [ $(diff $tmp1 $tmp2 | wc -l) -ne 0 ]
then
    echo "non-equal read names"
    exit 1
fi

# Saving the result.
cut -f 1-2,4 ${b1}.featureCounts > $tmp1
cut -f 2,4 ${b2}.featureCounts > $tmp2
paste $tmp1 $tmp2 | gzip > ${output}/${prefix}.raw.gz

## Doing some processing for Zika.
#zcat ${output}/${prefix}.raw.gz | grep "Zika" | grep "Assigned.*Assigned" | cut -f1,3,5 | gzip > ${output}/${prefix}.filtered.gz

rm -r ${prefix}
