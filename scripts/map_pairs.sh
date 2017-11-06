####################################
## This script extracts the UMIs, aligns each read to the hg38+Zika genome and merges the BAM files together.

set -e 
set -u

index=/lustre/jmlab/resources/genomes/star/hg38_zikv/
umilen=4
staropts=""
while getopts "i:o:u:" opt
do
    case $opt in
    i) # get new index
        index=$OPTARG
        ;;
    o) # get STAR options
        staropts=$OPTARG
        ;;
    u) # get UMI length
        umilen=$OPTARG
        ;;
    esac
done
shift "$((OPTIND - 1))"

if [ $# -ne 3 ]
then
    echo "$0 [-i <INDEX> -u <UMI LENGTH> -o <STAR OPTIONS>] <LIB1> <LIB2> <PREFIX>"
    exit 1
fi

if [ ! -e bam ]
then
    mkdir bam
    StripeThis bam
fi

if [ ! -e logs ]
then 
    mkdir logs
fi

first=$1
second=$2 
prefix=$3

if [ ! -e $prefix ]
then
    mkdir $prefix
    StripeThis $prefix
fi

####################################
# Running UMI-tools to extract the UMI from the first four bases of read 2.
new1=${prefix}/$(basename $first)
new2=${prefix}/$(basename $second)
bcpattern=$(printf '%*s' "${umilen}" | tr ' ' "N")
zcat $second | umi_tools extract --bc-pattern=${bcpattern} --read2-in=$first --read2-out=$new1 -L logs/${prefix}.umi-extract.log | gzip > $new2

####################################
# Running STAR to align to the genome.
# Done separately for each read to avoid favouring adjacent mapping locations for paired reads.

location=$( dirname $BASH_SOURCE )
for i in 1 2
do 
    if [ $i -eq 1 ]
    then 
        current=${new1}
        superfix=first
    else
        current=${new2}
	    superfix=second
    fi

    # Aligning:
    # --outFilterMultimapNmax specifies only unique alignments
    # --chimSegmentMin requires at least 20 bp to consider a chimeric alignment
    # --chimScoreJunctionNonGTAG removes the restriction for splice-site chimeras
    # --chimMainSegmentMultNmax only considers unique chimeric segments
    STAR --runThreadN 12 --genomeDir ${index} \
        --readFilesIn ${current} --readFilesCommand zcat \
        --outFileNamePrefix ${prefix}/${superfix}_ --outSAMtype BAM Unsorted --outSAMunmapped Within \
        --outFilterMultimapNmax 1 \
        --chimOutType WithinBAM \
        --chimSegmentMin 20 \
        --chimScoreJunctionNonGTAG 0 \
        --chimMainSegmentMultNmax 1 \
        ${staropts} 

    # Pulling out primary alignments and flipping the flag.
    curfile=${prefix}/${superfix}_Aligned.out.bam
    newfile=${prefix}/${superfix}.bam
    python ${location}/switch_flags.py $curfile $newfile ${superfix}
done

####################################
# Merging everyone together.

subfinal=${prefix}/out.bam
prefinal=${prefix}/whee
final=bam/${prefix}.bam
samtools cat -o ${subfinal} ${prefix}/first.bam ${prefix}/second.bam
samtools collate ${subfinal} ${prefinal}
samtools fixmate ${prefinal}.bam ${final}

rm -r ${prefix}
