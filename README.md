# Supporting scripts and data for "A root-specific NLR network confers resistance to plant parasitic nematodes"

Resources:

Software                            | Source
------------------------------------| ------------------------------------
*MAFFT v7.520*                      | (https://github.com/GSLBiotech/mafft)
*FastTree v2.1.11*                  | (http://www.microbesonline.org/fasttree/)
*Dendroscope v3.8.8*                | (https://software-ab.cs.uni-tuebingen.de/download/dendroscope3/welcome.html)
*R v4.3.1*                          | (https://cran.r-project.org/)
*NLRtracker*                        | (https://github.com/slt666666/NLRtracker)
*gene_cluster_matrix v0.1.4*        | (https://github.com/slt666666/gene-cluster-matrix)
*HISAT2 v2.1.0*                     | (https://daehwankimlab.github.io/hisat2/)
*STRINGTIE2 v2.1.1*                 | (https://github.com/gpertea/stringtie/releases)
*SAMtools v1.12*                    | (https://github.com/samtools/samtools/releases/)
*GFFcompare v0.11.6*                | (https://github.com/gpertea/gffcompare/releases)
*python v3.7.2*                     | (https://www.python.org/downloads/release/python-372/)

R packages:
```R
install.packages("tidyverse")
install.packages("readxl")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biostrings")
BiocManager::install("ggtree")
```
