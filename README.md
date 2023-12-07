# Supporting scripts and data for "A root-specific NLR network confers resistance to plant parasitic nematodes"

Resources:

Software                            | Source
------------------------------------| ------------------------------------
*MAFFT v7.520*                      | (https://github.com/GSLBiotech/mafft)
*FastTree v2.1.11*                  | (http://www.microbesonline.org/fasttree/)
*Dendroscope v3.8.8*                | (https://software-ab.cs.uni-tuebingen.de/download/dendroscope3/welcome.html)
*R v4.3.1*                          | (https://cran.r-project.org/)
*NLRtracker*                        | (https://github.com/slt666666/NLRtracker)

R packages:
```R
install.packages("tidyverse")
install.packages("readxl")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biostrings")
BiocManager::install("ggtree")
```
