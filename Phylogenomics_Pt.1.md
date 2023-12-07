Install and load the required packages:

```R
# Required packages
library(tidyverse)
library(Biostrings)
library(readxl)
library(ggtree)
```

Import the required reference files:

```R
# import genome data
Genome_data <- read_csv("/path/to/table_S2.csv")

## import the RefPlantNLR sequences
RefPlantNLR <- readAAStringSet("/path/to/RefPlantNLR.fasta")
RefPlantNLR_NBARC <- readAAStringSet("/path/to/RefPlantNLR_NBARC.fasta")

## import the Hero Cluster NLRs (HCN)
HCN <- readAAStringSet("/path/to/HCN.fasta")
HCN_NBARC <- readAAStringSet("/path/to/HCN_NBARC.fasta")

## import the tomato NRCs
SlNRC <- readAAStringSet("/path/to/Sl_NRC.fasta")
SlNRC_NBARC <- readAAStringSet("/path/to/Sl_NRC_NBARC.fasta")


## import the pepper NRCs
CaNRC <- readAAStringSet("/path/to/Ca_NRC.fasta")
CaNRC_NBARC <- readAAStringSet("/path/to/Ca_NRC_NBARC.fasta")
```

Now import the domain output from NLRtracker and extract the NLRs and other required data:
```R
## import the sequences
# Set the working directory to where you download the NLRtracker outputs
setwd("/path/to/NLRtracker")

# Make a function to add the file name to the proteins
read_delim_name <- function(flnm) {
  read_delim(flnm) %>%
    mutate(filename = flnm)
}


# Import the metadata
metadata <- list.files(pattern = "*_Domains.tsv", # metadata
                       full.names = T, recursive = TRUE) %>%
  map_df(~read_delim_name(.)) # bind all .tsv files together in a single dataframe
names(metadata)[11] <- "file_name" # correct the column name


# Modify the filename to RefSeq ID
metadata$file_name <- gsub("*_Domains.tsv", "", metadata$file_name)
metadata$file_name <- gsub(".*/", "", metadata$file_name)


# Subset based on "Status"
NLR_meta <- filter(metadata, metadata$Status == "NLR") # NLRs
NLR_meta <- filter(NLR_meta, NLR_meta$type == "CHAIN") # full length NLRs


# filter the unwanted architectures
NLR_meta_filtered <- NLR_meta[NLR_meta$Simple %in% c("CNL","BCNL"),]

# subset NBARC
NBARC_meta_filterd <- metadata[metadata$seqname %in% NLR_meta_filtered$seqname,]
NBARC_meta_filterd <- NBARC_meta_filterd[NBARC_meta_filterd$description == "NBARC",]


# convert NLR sequences dataframe to a Biostring object
NLR_meta_filtered_seq <- AAStringSet(NLR_meta_filtered$sequence)
NLR_meta_filtered_seq@ranges@NAMES <- NLR_meta_filtered$seqname

NBARC_meta_filtered_seq <- AAStringSet(NBARC_meta_filterd$sequence)
NBARC_meta_filtered_seq@ranges@NAMES <- NBARC_meta_filterd$seqname


writeXStringSet(NBARC_meta_filtered_seq, "/path/to/NBARC_seq.csv")
writeXStringSet(NLR_meta_filtered_seq, "/path/to/NLR_seq.csv")


# merge the sequences with reference sequences
NLR_ref_seq <- c(NLR_meta_filtered_seq, RefPlantNLR, HCN, SlNRC, CaNRC)
NBARC_ref_seq <- c(NBARC_meta_filtered_seq, RefPlantNLR_NBARC, HCN_NBARC, SlNRC_NBARC, CaNRC_NBARC)
```

Filter out truncated or unusually long NBARC sequences. We later use this for building phylogenetic trees:
```R
# filter out truncated NBARC sequences

NBARC_ref_seq_filtered <- NBARC_ref_seq[NBARC_ref_seq@ranges@width > 300] # remove anything shorter than 300 amino acids
NBARC_ref_seq_filtered <- NBARC_ref_seq_filtered[NBARC_ref_seq_filtered@ranges@width < 400] # remove anything longer than 400 amino acids


# export the NBARC file to align and make a phylogenetic tree

writeXStringSet(NBARC_ref_seq, "/path/to/NBARC_ref_hero.fasta")
writeXStringSet(NBARC_ref_seq_filtered, "/path/to/NBARC_ref_hero_filtered.fasta")

Using the exported NBARC sequences we make a tree with the following two steps:
```bash
# alignment
mafft --anysymbol NBARC_ref_hero.fasta > NBARC_ref_hero.fasta.afa

# tree construction
FastTree -lg NBARC_ref_hero.afa > NBARC_ref_hero_lg.newick
```



