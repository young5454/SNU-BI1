# ------------------------------------------------------------------------------
# Script to conduct GO analysis of CLIP-seq & Ribo-footprinting data [052824]
# ------------------------------------------------------------------------------
# 1. Load libraries & Set workspace
# ------------------------------------------------------------------------------
library(biomaRt)
library(latex2exp)
library(org.Mm.eg.db)
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

filt.meta.counts <- meta.counts[rna.counts >= 30 & ribo.counts >= 80, ]
nrow(filt.meta.counts)  # 55158 -> 8154

# ------------------------------------------------------------------------------
# 4. Calculate CLIP enrichment & Rden change
# ------------------------------------------------------------------------------
## CLIP-enrichment
### Normalize by dividing with total read sum of filtered counts
norm.clip.35L33G <- filt.meta.counts$CLIP.35L33G / sum(filt.meta.counts$CLIP.35L33G)
norm.rna.control <- filt.meta.counts$RNA.control / sum(filt.meta.counts$RNA.control)
clip.enrichment <- log2(norm.clip.35L33G / norm.rna.control)  

## Rden change
### Normalize by dividing with total read sum of filtered counts
density.silin28a <- filt.meta.counts$RPF.siLin28a / filt.meta.counts$RNA.siLin28a
density.siluc <- filt.meta.counts$RPF.siLuc / filt.meta.counts$RNA.siLuc
rden.change <- log2(density.silin28a / density.siluc)

## Merge data to counts
filt.meta.counts$CLIP.Enrichment <- clip.enrichment
filt.meta.counts$Rden.Change <- rden.change
clip.rden2 <- data.frame(clip.enrichment, rden.change)

## Save filtered counts
write.table(filt.meta.counts, "./SNU-BI1.Data/filt.meta.counts.csv", sep=",")

## Plot - Fig4D
scatter <- ggplot(clip.rden2, aes(x=clip.enrichment, y=rden.change)) + 
            geom_point(size=0.2, alpha=0.25) + 
            labs(title=TeX("CLIP and ribosome footrprinting upon $\\textit{Lin28a}$ knockdown"),
                 x=TeX("LIN28A CLIP enrichment ($\\log_2$)"), 
                 y=TeX("Ribosome density change upon $\\textit{Lin28a}$ knockdown ($\\log_2$)")) +
            coord_cartesian(xlim=c(-6, 4), ylim=c(-3, 2), expand=FALSE) +
            theme_classic() +
            theme(
              panel.grid.major=element_line(color="gray", size=0.25),
              panel.background=element_rect(fill="white"),
              plot.background=element_rect(fill="white"),
              plot.title=element_text(family="Arial", size=9.5),
              axis.title.x=element_text(family="Arial", size=9),
              axis.title.y=element_text(family="Arial", size=9),
              axis.text.x=element_text(size=11),
              axis.text.y=element_text(size=11),
              axis.ticks=element_line(size=0.25),
              axis.line=element_line(size=0.25),
              axis.line.x.top=element_blank(),      # Remove top axis line
              axis.line.y.right=element_blank())
ggsave(scatter, filename="./SNU-BI1/FreeExpansionPlots/fig4d.png", 
       width=4, height=4, units='in', dpi=600)

# ------------------------------------------------------------------------------
# 4. Map GO terms
# ------------------------------------------------------------------------------
## Define gene list
filt.gene.list <- filt.meta.counts$GeneID

## Download Mapping Table - 2024.05.29 last mod.
## Download GOA DB - 2024.04.16 last mod.
mus.goa <- readLines("./SNU-BI1.Data/goa_mouse.gaf")
mus.goa <- mus.goa[!grepl("^!", mus.goa)]
mus.goa <- strsplit(mus.goa, "\t")
mus.goa.df <- do.call(rbind, mus.goa)
mus.goa.df <- as.data.frame(mus.goa.df)
mus.goa.df <- subset(mus.goa.df, select=-ncol(mus.goa.df))

## Rename GOA columns
colnames(mus.goa.df) <- c("DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier",
                          "GO_ID", "DB:Reference", "Evidence_Code", "With_From",
                          "Aspect", "DB_Object_Name", "DB_Object_Synonym", 
                          "DB_Object_Type", "Taxon", "Date", "Assigned_By")

## Map MGIsymbol to GO term
go.matched <- mus.goa.df[mus.goa.df$DB_Object_Symbol %in% filt.meta.counts$MGIsymbol, ]
mapping.table <- go.matched[, c(3, 5)]   # Total 7706 / 8154 genes are mapped






