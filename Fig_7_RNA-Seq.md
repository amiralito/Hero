# RNA-Seq pipeline - read alignment to genome

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
# RNA-Seq pipeline - quantification

Assemble and quantify expressed genes and transcripts:

```bash
stringtie -p 8 -G annotation.gtf -o outputfile.gtf -l filename /path/to/bam/files/file.sorted.bam
```

Merge for all samples (mergelist_stringtie.txt needs to contain the name of the outputfile.gtf files):

```bash
stringtie --merge -p 8 -G annotation.gtf -o output_stringtie.gtf mergelist_stringtie.txt
```

Compare gene annotations:

```bash
gffcompare -r annotation.gtf -G -o merged output_stringtie.gtf

stringtie -e -B -p 8 -G output_stringtie.gtf -o /path/to/merge/file.gtf /path/to/bam/files/file.sorted.bam
```

Prepare count matrix for Deseq2 (sample_list.txt needs to contain paths and names of .gtf files):

```bash
python prepDE3.py -i sample_lst.txt
```
