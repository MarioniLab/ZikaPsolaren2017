# Mapping libraries
mapscript=../scripts/map_pairs.sh
bsub -R "rusage[mem=40000]" -n 12 -e log1.err -o log1.err bash ${mapscript} fastq/Lib2_S1_L001_R1_001.fastq.gz fastq/Lib2_S1_L001_R2_001.fastq.gz Lib2_L001
bsub -R "rusage[mem=40000]" -n 12 -e log2.err -o log2.err bash ${mapscript} fastq/Lib2_S1_L002_R1_001.fastq.gz fastq/Lib2_S1_L002_R2_001.fastq.gz Lib2_L002
bsub -R "rusage[mem=40000]" -n 12 -e log1c.err -o log1c.err bash ${mapscript} fastq/Lib2-cont_S2_L001_R1_001.fastq.gz fastq/Lib2-cont_S2_L001_R2_001.fastq.gz Lib2-cont_L001
bsub -R "rusage[mem=40000]" -n 12 -e log2c.err -o log2c.err bash ${mapscript} fastq/Lib2-cont_S2_L002_R1_001.fastq.gz fastq/Lib2-cont_S2_L002_R2_001.fastq.gz Lib2-cont_L002

# Dedupping libraries
dupscript=../scripts/dedup_pairs.sh
bsub -R "rusage[mem=16000]" -n 1 -e logd.err -o logd.err bash ${dupscript} Lib2 bam/Lib2_L001.bam bam/Lib2_L002.bam
bsub -R "rusage[mem=10000]" -n 1 -e logdc.err -o logdc.err bash ${dupscript} Lib2-cont bam/Lib2-cont_L001.bam bam/Lib2-cont_L002.bam

# Building HDF5 libraries.
buildscript=../scripts/build_hdf5.sh
bsub -R "rusage[mem=10000]" -n 1 -e log1.err -o log1.err bash ${buildscript} bam/Lib2.bam processed/Lib2.h5
bsub -R "rusage[mem=10000]" -n 1 -e log2.err -o log2.err bash ${buildscript} bam/Lib2-cont.bam processed/Lib2-cont.h5

