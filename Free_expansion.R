# ------------------------------------------------------------------------------
# Script to conduct GO analysis of CLIP-seq & Ribo-footprinting data [052824]
# ------------------------------------------------------------------------------
# 1. Load libraries & Set workspace
# ------------------------------------------------------------------------------
library(biomaRt)
library(clusterProfiler)
library(ggplot2)
library(colorRamp2)

setwd("/home/local/hoeyoungkim_000504/coursework/")

# ------------------------------------------------------------------------------
# 2. Match annotation to raw counts
# ------------------------------------------------------------------------------
## Load raw counts data generated with FeatureCounts
raw.counts <- read.csv("./SNU-BI1.Data/raw.counts.txt", sep="\t", header=TRUE)

## Rename header row
header <- c("GeneID", "Chr", "Start", "End", "Strand", "Length",
            "CLIP.35L33G", "RNA.control", 
            "RNA.siLin28a", "RNA.siLuc",
            "RPF.siLin28a", "RPF.siLuc")
colnames(raw.counts) <- header

## Remove Chr, Start, End, Strand info
raw.counts <- raw.counts[, -c(2, 3, 4, 5)]
head(raw.counts)

## Check for duplicates - 55359 unique rows
gene_ids <- raw.counts$GeneID
duplicates <- gene_ids[duplicated(gene_ids)]
duplicates   # 0

gene_ids <- sub("\\..*", "", gene_ids)
duplicates <- gene_ids[duplicated(gene_ids)]
duplicates   # 0

raw.counts$GeneID <- gene_ids

## Match MGI-symbol and Description using BiomaRt
mus.ensembl <- useEnsembl("ensembl", dataset="mmusculus_gene_ensembl", version=111)
mus.meta <- getBM(attributes=c("ensembl_gene_id", "mgi_symbol", "description"),
                  filters="ensembl_gene_id",
                  values=gene_ids,
                  mart=mus.ensembl)
colnames(mus.meta) <- c("GeneID", "MGIsymbol", "Description")

## Add matched MGI-symbol column 
matched.counts <- raw.counts[raw.counts$GeneID %in% mus.meta$GeneID, ]
order <- match(matched.counts$GeneID, mus.meta$GeneID)
mus.meta <- mus.meta[order, ]
meta.counts <- cbind(matched.counts[, c(1,2)],
                     mus.meta[, c(2,3)],
                     matched.counts[, 2:ncol(meta.filter.counts)])
## Save Meta-Counts
write.table(meta.counts, "./SNU-BI1.Data/meta.counts.csv", sep=",")

# ------------------------------------------------------------------------------
# 3. Filter counts 
# ------------------------------------------------------------------------------
## Filter 1. <30 reads in RNA-seq
rna.counts <- rowSums(meta.counts[, c(7,8,9)])

## Filter 2. <80 raw footprint tags in siLuc library
ribo.counts <- meta.counts$RPF.siLuc

filtered.meta.counts <- meta.counts[rna.counts >= 30 & ribo.counts >= 80, ]
nrow(filtered.meta.counts)  # 55158 -> 8154
