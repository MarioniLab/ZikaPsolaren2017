# Mapping libraries
mapscript=../scripts/map_pairs.sh
for i in {1..18}
do
    stub=TestLib8-${i}_S${i}
    bsub -R "rusage[mem=40000]" -n 12 -e log${i}-1.err -o log${i}-1.out bash ${mapscript} fastq/${stub}_L001_R1_001.fastq.gz fastq/${stub}_L001_R2_001.fastq.gz ${stub}_L001
    bsub -R "rusage[mem=40000]" -n 12 -e log${i}-2.err -o log${i}-2.out bash ${mapscript} fastq/${stub}_L002_R1_001.fastq.gz fastq/${stub}_L002_R2_001.fastq.gz ${stub}_L002
done

# Dedupping libraries
dupscript=../scripts/dedup_pairs.sh
for i in {1..18}
do
    stub=TestLib8-${i}_S${i}
    bsub -R "rusage[mem=10000]" -n 1 -e log${i}.err -o log${i}.out bash ${dupscript} ${stub} bam/${stub}_L001.bam bam/${stub}_L002.bam
done

# Counting libraries.
countscript=../scripts/count_pairs.sh
for i in {1..18}
do
    stub=TestLib8-${i}_S${i}
    bsub -R "rusage[mem=10000]" -n 1 -e log${i}.err -o log${i}.out bash ${countscript} bam/${stub}.bam ../scripts/combined.gtf processed/
done

# Building HDF5 libraries.
buildscript=../scripts/build_hdf5.sh
for i in {1..18}
do
    stub=TestLib8-${i}_S${i}
    bsub -R "rusage[mem=10000]" -n 1 -e log${i}.err -o log${i}.out bash ${buildscript} bam/${stub}.bam processed/${stub}.h5
done
