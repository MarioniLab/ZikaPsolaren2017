## This script merges technical replicates and removes UMI-labelled duplicates.

set -e 
set -u

if [ $# -le 1 ]
then        
    echo "$0 <PREFIX> <INPUT FILES...>"
    exit 1
fi

if [ ! -e bam ]
then
    mkdir bam
    StripeThis bam
fi

prefix=$1

if [ ! -e $prefix ]
then
    mkdir $prefix
    StripeThis $prefix
fi

shift
NEWFILES=()
for f in "$@"
do
    newfile=${prefix}/$(basename ${f})
    NEWFILES+=($newfile)
    tmp=${prefix}/temp.bam
    samtools view -h -F 2048 -b ${f} -o ${tmp}
    samtools sort -o $newfile ${tmp}
    rm ${tmp}
done

semifinal=${prefix}/out.bam
samtools merge ${semifinal} "${NEWFILES[@]}"
samtools index ${semifinal}

prefinal=${prefix}/dedup.bam
umi_tools dedup -I ${semifinal} --paired -S ${prefinal}
#MarkDuplicates I=${semifinal} O=${prefinal} M=${prefix}/blah.log REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT

final=bam/${prefix}.bam
samtools sort -n -o ${final} ${prefinal}

rm -r ${prefix}
