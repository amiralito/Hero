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

# RNA-Seq pipeline - data analysis using DeSeq2

Load required libraries:
```R
library("DESeq2")
library("dplyr")
library("ggplot2")
library("reshape2")
library("e1071")
```

Input data from DeSeq2:
```R
# Load raw data
countData <- as.matrix(read.csv("/path/to/folder/gene_count_matrix.csv", row.names="gene_id"))
colData <- read.csv("/path/to/pheno/data/PHENO_DATA.csv", row.names=1)

# Check that all sample IDs in colData are also in countData and match their orders
# It's important to ensure data consistency between countData and colData
all(rownames(colData) %in% colnames(countData))
countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))

# Create a DESeqDataSet from count matrix and labels
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ tissue)
```

Run the default analysis for DESeq2"
```R
# This step estimates size factors, dispersions, and fits the negative binomial GLM
dds <- DESeq(dds)

# Extract results from the DESeq analysis for the entire dataset
res <- results(dds)
summary(res)

#Sort by adjusted p-value and display
(resOrdered <- res[order(res$padj), ])

#Extract number of differentially expressed genes with adjusted P < 0.05
res05 <- results(dds, alpha=0.05)
summary(res05)

# Extract normalized counts from the DESeqDataSet
table_counts_normalized <- counts(dds, normalized=TRUE)

# Write the normalized counts to a CSV file
write.csv(table_counts_normalized, file = "/path/to/folder/gene_count_matrix_normalized.csv", row.names = TRUE)

# Write the DESeq results to a CSV file
write.csv(res, file = "/path/to/folder/gene_count_matrix_DeSeq2.csv", row.names = TRUE)
write.csv(res05, file = "/path/to/folder/gene_count_matrix_DeSeq2_05.csv", row.names = TRUE)

#subsetting with partial ID match from provided .csv list with no header

# Read in and assign to variable for All_NLRs
All_NLRs <- read.csv("/path/to/folder/All_NLRs.csv", header = FALSE)
ALL_NLRs_sub <- c(All_NLRs$V1)
All_NLRs_sub_res <- subset(res, grepl(paste(ALL_NLRs_sub, collapse= "|"), rownames(res)))

# Write the subsetted results to a CSV file for All_NLRs
write.csv(All_NLRs_sub_res, file = "/path/to/folder/ALL_NLRs_res.csv", row.names = TRUE)
```
Generate figures for log2FoldChange data:
```R
# Create a histogram for the log2FoldChange values from the entire dataset
hist(res$log2FoldChange, xlim=c(-22,22), lwd=4, breaks=1000, lend=100, col="red", border="red")
# Add a vertical line at log2FoldChange = 0 for reference
abline(v = 0, col="white", lwd=1, lty=2)

# Create a histogram for the log2FoldChange values from the All_NLRs subset
hist(All_NLRs_sub_res$log2FoldChange, xlim=c(-22,22), lwd=2, breaks=20, lend=100, col="red", border="red")
# Add a vertical line at log2FoldChange = 0 for reference
abline(v = 0, col="white", lwd=1, lty=2)

# Calculate mean and skewness based on log2FoldChange
# Extract the log2FoldChange column
log2_fold_change <- All_NLRs_sub_res$log2FoldChange

# Check for "NA" values in log2FoldChange
na_values <- sum(is.na(log2_fold_change))

if (na_values > 0) {
  # If there are "NA" values, remove them before calculating skewness
  log2_fold_change <- log2_fold_change[!is.na(log2_fold_change)]
  
  # Calculate the mean
  mean_value <- mean(log2_fold_change)
  
  # Print the result
  cat("Mean:", mean_value, "\n")
  
  # Calculate skewness after removing "NA" values
  skewness_value <- skewness(log2_fold_change)
  
  # Print the skewness value
  print(skewness_value)
```
Generate figures for normalized counts data:
```R
# Subset the normalized counts for the All_NLRs subset
# Write the subsetted counts to a CSV file for the All_NLRs subset
All_NLRs_sub_counts <- subset(table_counts_normalized, grepl(paste(ALL_NLRs_sub, collapse= "|"), rownames(table_counts_normalized)))
write.csv(All_NLRs_sub_counts, file = "/path/to/folder/ALL_NLRs_counts.csv", row.names = TRUE)

#single gene plot counts examples
#HCN_cluster
plotCounts(dds, gene = "MSTRG.7595|Hero_A", intgroup = "tissue", returnData = TRUE) %>%
  ggplot() +
  aes(tissue, count, color = tissue) +
  geom_jitter(width = 0.2, size = 4, shape = 1, stroke = 3) +
  ggtitle("Hero_A") +
  theme_bw() +
  scale_y_log10(limits = c(1, 10000), oob = scales::squish)

#NRCs
plotCounts(dds, gene = "MSTRG.5251|Solyc03g005660.3", intgroup = "tissue", returnData = TRUE) %>%
  ggplot() +
  aes(tissue, count, color = tissue) +
  geom_jitter(width = 0.08, size = 4, shape = 1, stroke = 3) +
  ggtitle("NRCX") +
  theme_bw() +
  scale_y_log10(limits = c(1, 91000), oob = scales::squish)
```
