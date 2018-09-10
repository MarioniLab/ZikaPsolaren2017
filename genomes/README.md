# Generating annotation

Annotation for the hg38 build of the human genome was obtained from Ensembl and processed as shown below.

```sh
zcat Homo_sapiens.GRCh38.90.gtf.gz | \
    sed -r "s/^([0-9MXY])/chr\1/" | \
    sed "s/^chrMT/chrM/g" | \
    awk '$3 == "exon"' > hg38.gtf
```

The Zika virus GTF consists of a single entry spanning the start and end of the Zika genome, required for _STAR_ to treat it as a transcript.

# Generating STAR indices

Builds using the standard hg38 sequence and a Zika genome sequence obtained from Omer Ziv.

```{r}
mkdir temp
cat hg38.gtf ZIKV_H.sapiens_Brazil_PE243_2015-1.gtf > combined.gtf
mkdir hg38_zikv
bsub -R "rusage[mem=40000]" -n 4 -o log.out -e log.err STAR --runMode genomeGenerate --runThreadN 4 --genomeDir hg38_zikv \
    --genomeFastaFiles hsa.hg38.fa ZIKV_H.sapiens_Brazil_PE243_2015-1.fa \
    --sjdbGTFfile combined.gtf --sjdbOverhang 100
```

**Unused:** builds using non-coding RNAs (need `--genomeSAindexNbases` to avoid index bug).
Uses a manually curated set of human non-coding RNAs.

```{r}
mkdir hs_ncRNAs
bsub -R "rusage[mem=40000]" -n 4 -o log.out -e log.err STAR --runMode genomeGenerate --runThreadN 4 --genomeDir hs_ncRNAs \
     --genomeFastaFiles hs_ncRNAs.fa --genomeSAindexNbases 5
```

