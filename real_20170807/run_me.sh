# Mapping libraries (UMI of length 6)
mapscript=../scripts/map_pairs.sh
for i in $(ls fastq/*_R1_*)
do
    prefix=$(basename $i | sed "s/_R1.*//")
    other=$(echo $i | sed "s/R1/R2/")
    bsub -R "rusage[mem=40000]" -n 12 -e log-${prefix}.err -o log-${prefix}.out bash ${mapscript} -u 6 ${i} ${other} ${prefix} 
done

# Dedupping libraries
dupscript=../scripts/dedup_pairs.sh
for i in $(ls bam/*_L001.bam)
do
    prefix=$(basename $i | sed "s/_L001.*//")
    other=$(echo $i | sed "s/_L001/_L002/")
    bsub -R "rusage[mem=10000]" -n 1 -e log-${prefix}.err -o log-${prefix}.out bash ${dupscript} ${prefix} ${i} ${other}
done

# Counting libraries.
countscript=../scripts/count_pairs.sh
for i in $(ls bam/*.bam | grep "S[0-9].bam")
do
    prefix=$(basename $i | sed "s/.bam$//")
    bsub -R "rusage[mem=10000]" -n 1 -e log-${prefix}.err -o log-${prefix}.out bash ${countscript} ${i} ../genomes/combined.gtf processed/
done

# Building HDF5 libraries.
buildscript=../scripts/build_hdf5.sh
for i in $(ls bam/*.bam | grep "S[0-9].bam")
do
    prefix=$(basename $i | sed "s/.bam$//")
    bsub -R "rusage[mem=10000]" -n 1 -e log-${prefix}.err -o log-${prefix}.out bash ${buildscript} bam/${prefix}.bam processed/${prefix}.h5
done
