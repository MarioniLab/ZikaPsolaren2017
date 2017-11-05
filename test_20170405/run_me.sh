# Mapping libraries
mapscript=../scripts/map_pairs.sh
bsub -R "rusage[mem=40000]" -n 12 -e log1.err -o log1.err bash ${mapscript} fastq/Lib3-1_S1_L001_R1_001.fastq.gz fastq/Lib3-1_S1_L001_R2_001.fastq.gz Lib3-1_S1
bsub -R "rusage[mem=40000]" -n 12 -e log1.err -o log1.err bash ${mapscript} fastq/Lib3-1_S1_L002_R1_001.fastq.gz fastq/Lib3-1_S1_L002_R2_001.fastq.gz Lib3-1_S2

bsub -R "rusage[mem=40000]" -n 12 -e log2.err -o log2.err bash ${mapscript} fastq/Lib3-2_S2_L001_R1_001.fastq.gz fastq/Lib3-2_S2_L001_R2_001.fastq.gz Lib3-2_S1  
bsub -R "rusage[mem=40000]" -n 12 -e log2.err -o log2.err bash ${mapscript} fastq/Lib3-2_S2_L002_R1_001.fastq.gz fastq/Lib3-2_S2_L002_R2_001.fastq.gz Lib3-2_S2

bsub -R "rusage[mem=40000]" -n 12 -e log3.err -o log3.err bash ${mapscript} fastq/Lib3-3_S3_L001_R1_001.fastq.gz fastq/Lib3-3_S3_L001_R2_001.fastq.gz Lib3-3_S1
bsub -R "rusage[mem=40000]" -n 12 -e log3.err -o log3.err bash ${mapscript} fastq/Lib3-3_S3_L002_R1_001.fastq.gz fastq/Lib3-3_S3_L002_R2_001.fastq.gz Lib3-3_S2

bsub -R "rusage[mem=40000]" -n 12 -e log4.err -o log4.err bash ${mapscript} fastq/Lib3-4_S4_L001_R1_001.fastq.gz fastq/Lib3-4_S4_L001_R2_001.fastq.gz Lib3-4_S1
bsub -R "rusage[mem=40000]" -n 12 -e log4.err -o log4.err bash ${mapscript} fastq/Lib3-4_S4_L002_R1_001.fastq.gz fastq/Lib3-4_S4_L002_R2_001.fastq.gz Lib3-4_S2

bsub -R "rusage[mem=40000]" -n 12 -e log6.err -o log6.err bash ${mapscript} fastq/Lib3-6_S5_L001_R1_001.fastq.gz fastq/Lib3-6_S5_L001_R2_001.fastq.gz Lib3-6_S1
bsub -R "rusage[mem=40000]" -n 12 -e log6.err -o log6.err bash ${mapscript} fastq/Lib3-6_S5_L002_R1_001.fastq.gz fastq/Lib3-6_S5_L002_R2_001.fastq.gz Lib3-6_S2

# Dedupping libraries
dupscript=../scripts/dedup_pairs.sh
bsub -R "rusage[mem=10000]" -n 1 -e log1.err -o log1.err bash ${dupscript} Lib3-1 bam/Lib3-1_S1.bam bam/Lib3-1_S2.bam
bsub -R "rusage[mem=10000]" -n 1 -e log2.err -o log2.err bash ${dupscript} Lib3-2 bam/Lib3-2_S1.bam bam/Lib3-2_S2.bam
bsub -R "rusage[mem=16000]" -n 1 -e log3.err -o log3.err bash ${dupscript} Lib3-3 bam/Lib3-3_S1.bam bam/Lib3-3_S2.bam
bsub -R "rusage[mem=10000]" -n 1 -e log4.err -o log4.err bash ${dupscript} Lib3-4 bam/Lib3-4_S1.bam bam/Lib3-4_S2.bam
bsub -R "rusage[mem=10000]" -n 1 -e log6.err -o log6.err bash ${dupscript} Lib3-6 bam/Lib3-6_S1.bam bam/Lib3-6_S2.bam

# Counting libraries
countscript=../scripts/count_pairs.sh
bsub -R "rusage[mem=5000]" -n 1 -o log1.out -e log1.err bash ${countscript} bam/Lib3-1.bam ../scripts/combined.gtf processed
bsub -R "rusage[mem=5000]" -n 1 -o log2.out -e log2.err bash ${countscript} bam/Lib3-2.bam ../scripts/combined.gtf processed
bsub -R "rusage[mem=5000]" -n 1 -o log3.out -e log3.err bash ${countscript} bam/Lib3-3.bam ../scripts/combined.gtf processed
bsub -R "rusage[mem=5000]" -n 1 -o log4.out -e log4.err bash ${countscript} bam/Lib3-4.bam ../scripts/combined.gtf processed
bsub -R "rusage[mem=5000]" -n 1 -o log6.out -e log6.err bash ${countscript} bam/Lib3-6.bam ../scripts/combined.gtf processed

# Building HDF5 libraries.
buildscript=../scripts/build_hdf5.sh
bsub -R "rusage[mem=10000]" -n 1 -e log1.err -o log1.err bash ${buildscript} bam/Lib3-1.bam processed/Lib3-1.h5
bsub -R "rusage[mem=10000]" -n 1 -e log2.err -o log2.err bash ${buildscript} bam/Lib3-2.bam processed/Lib3-2.h5
bsub -R "rusage[mem=10000]" -n 1 -e log3.err -o log3.err bash ${buildscript} bam/Lib3-3.bam processed/Lib3-3.h5
bsub -R "rusage[mem=10000]" -n 1 -e log4.err -o log4.err bash ${buildscript} bam/Lib3-4.bam processed/Lib3-4.h5
bsub -R "rusage[mem=10000]" -n 1 -e log6.err -o log6.err bash ${buildscript} bam/Lib3-6.bam processed/Lib3-6.h5

