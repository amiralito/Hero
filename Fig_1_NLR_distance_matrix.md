NLR annotations from the ITAG3.2 and RefSeq SL3.0 gff3 files using NLRtracker
```bash
# NLR annotation
# ITAG3.2_proteins.fasta ... https://solgenomics.net/ftp/genomes/Solanum_lycopersicum/annotation/ITAG3.2_release/ITAG3.2_proteins.fasta
# RefSeq.SL3.0.protein.faa ... download from https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000188115.4/
NLRtracker.sh -s ITAG3.2_proteins.fasta -o Tomato_NLRs_1
NLRtracker.sh -s RefSeq.SL3.0.protein.faa -o Tomato_NLRs_2
```

Using NBARC sequences to construct phylogenetic tree
```bash
# alignment
# Tomato_ITAG3.2_with_RefSeq.fasta ... merged amino acid sequences of NLRs that were annotated from ITAG3.2 & RefSeq SL3.0 and AtZAR1
mafft --anysymbol Tomato_ITAG3.2_with_RefSeq.fasta > Tomato_ITAG3.2_with_RefSeq_alignment.fasta

# After alignment
# extract NBARC sequences from Tomato_ITAG3.2_with_RefSeq_alignment.fasta based on AtZAR1
# NBARC sequences were manually edited to remove NLRs with truncated p-loop domain -> Tomato_ITAG3.2_with_RefSeq_NBARC_alignment.fasta

# tree construction by RAxML
raxmlHPC-PTHREADS-AVX2 -# 100 -f a -x 1024 -m PROTGAMMAAUTO -s Tomato_ITAG3.2_with_RefSeq_NBARC_alignment.fasta -p 121 -n Tomato_ITAG3.2_with_RefSeq
```

Construct distance matrix by gene_cluster_matrix
```bash
# NLR_positions.csv ... positions of NLRs from gff files; GCF_000188115.5_SL3.1_genomic.gff & ITAG3.2_gene_models.gff
# NLR_clade.csv ... clade information defined based on phylogenetic tree
# gene_ids.txt ... the order of gene ids based on phylogenetic tree 
gene_cluster_matrix -p NLR_positions.csv -f mRNA -c NLR_clade.csv -i gene_ids.txt -o Tomato_matrix
```
