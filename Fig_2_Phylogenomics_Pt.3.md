The phylogenetic trees extracted in each step in the [previous part](Fig_2_Phylogenomics_Pt.2.md) are imported in R and are used to extract the sequences belonging to each clade. The sequences are then exported to be realigned and generate the next phylogenetic trees:

NRC Superclade:
```R
# import the NRC superclade extracted based on phylogeny
NRC_tree <- read.tree("/path/to/NRC_superclade.tree")

# extract the NRC metadata and sequences
NRC_NBARC_seq <- NBARC_ref_seq_filtered[NBARC_ref_seq_filtered@ranges@NAMES %in% NRC_tree$tip.label,]

# export for realignment and making the tree
writeXStringSet(NRC_NBARC_seq, "/path/to/NRC_superclade_NBARC.fasta")
```

NRC Helpers:
```R
# import NRC Helper tree
NRCH_tree <- read.tree("/path/to/NRCH.tree")

# extract the NRC metadata and sequences
NRCH_NBARC_seq <- NBARC_ref_seq[NBARC_ref_seq@ranges@NAMES %in% NRCH_tree$tip.label,]

# export for realignment and making the tree
writeXStringSet(NRCH_NBARC_seq, "/path/to/NRCH_NBARC.fasta")
```

NRC SD-type sensors:
```R
# import NRC SD-type sensor tree
NRC_SD_tree <- read.tree("/path/to/NRC_SD.tree")

# extract the NRC metadata and sequences
NRC_SD_NBARC_seq <- NBARC_ref_seq[NBARC_ref_seq@ranges@NAMES %in% NRC_SD_tree$tip.label,]

# export for realignment and making the tree
writeXStringSet(NRC_SD_NBARC_seq, "/path/to/NRC_SD_NBARC.fasta")
```

Hero Clade:
```R
# import the Hero subclade from SD-type sensors
NRC_Hero_tree <- read.tree("/path/to/NRC_SD_Hero.tree")

# extract the NRC metadata and sequences
NRC_Hero_NBARC_seq <- NBARC_ref_seq[NBARC_ref_seq@ranges@NAMES %in% NRC_Hero_tree$tip.label,]

# export for realignment and making the tree
writeXStringSet(NRC_Hero_NBARC_seq, "/path/to/NRC_SD_Hero_NBARC.fasta")
```

NRC6 Clade:
```R
# import the NRC6 clade
NRC6_tree <- read.tree("~/Desktop/NRC6_hero/analysis/phylogeny/NRC6.tree")
```
No new trees were made for this clade.





Finally, calculate the frequency of NRC helpers, NRC6 clade, NRC SD-type sensors, and Hero Clade in the genomes used in this study:
```R
# extract metadata and calculate the frequency of each category based on species

NRC_SD_meta <- NLR_meta_filtered[NLR_meta_filtered$seqname %in% NRC_SD_tree$tip.label,]
NRC_SD_speceis_freq <- table(NRC_SD_meta$file_name) %>% as.data.frame() %>% setNames(c("file_name","SD"))


NRC_Hero_meta <- NLR_meta_filtered[NLR_meta_filtered$seqname %in% NRC_Hero_tree$tip.label,]
NRC_Hero_species_freq <- table(NRC_Hero_meta$file_name) %>% as.data.frame() %>% setNames(c("file_name","Hero"))


NRCH_meta <-  NLR_meta_filtered[NLR_meta_filtered$seqname %in% NRCH_tree$tip.label,]
NRCH_species_freq <- table(NRCH_meta$file_name) %>% as.data.frame() %>% setNames(c("file_name","NRCH"))


NRC6_meta <- NLR_meta_filtered[NLR_meta_filtered$seqname %in% NRC6_tree$tip.label,]
NRC6_species_freq <- table(NRC6_meta$file_name) %>% as.data.frame() %>% setNames(c("file_name","NRC6"))



# species dist dataframe

species_freq <- Genome_data %>% left_join(NRC_SD_speceis_freq, by = "file_name")
species_freq <- species_freq %>% left_join(NRC_Hero_species_freq, by = "file_name")
species_freq <- species_freq %>% left_join(NRCH_species_freq, by = "file_name")
species_freq <- species_freq %>% left_join(NRC6_species_freq, by = "file_name")


species_freq_filtered <- species_freq[,c(1,2,3,17,18,19,20)]

species_freq_filtered[is.na(species_freq_filtered)] <- 0

write_csv(species_freq_filtered, "/path/to/species_freq.csv")
```
