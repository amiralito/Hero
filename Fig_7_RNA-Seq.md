# RNA-Seq pipeline

Create transcriptome index for HISAT2:

```bash
hisat2-build -p 16 /path/to/genome/assembly/file.fasta /output/path/to/index/file
```

Use HISAT2 to align reads pairs:

```bash
hisat2 -q -p 16 --pen-canintronlen G,-8,1 --pen-noncanintronlen G,-8,1 --no-mixed --no-discordant --dta -x path/to/index/file -1 /path/to/first/readpair/file.fq -2 /path/to/second/readpair/file.fq -S output.sam
```

Use samtools to convert .sam to .bam, sort, and index:

```bash
samtools view -b file.sam > file.bam
samtools sort file.bam > file.sorted.bam
samtools index file.sorted.bam
```

