set -e 
set -u

if [ $# -ne 2 ]
then
    echo "$0 <BAM> <HDF5>"
    exit 1
fi

curfile=$1
output=$2
prefix=$(basename $curfile | sed "s/.bam$//")
if [ ! -e $prefix ]
then
    mkdir ${prefix}
    StripeThis ${prefix}
fi

# Stripping out supplementary alignments.
tmp=${prefix}/tmp.bam
samtools view -h -F 2048 -b ${curfile} -o ${tmp}

# Running through diffHic.
echo "
library(diffHic)
library(Rsamtools)
chrs <- scanBamHeader('${curfile}')[[1]][[1]]
ranges <- emptyGenome(chrs)
param <- pairParam(ranges)
out <- preparePairs('${curfile}', param, file='${output}', storage=1000, minq=255, dedup=FALSE)
save(param, out, file=sub('h5$', 'Rda', '${output}'))
" | Rdevel --no-save --vanilla

rm -r ${prefix}
