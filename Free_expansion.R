# ------------------------------------------------------------------------------
# Script to conduct GO analysis of CLIP-seq & Ribo-footprinting data [052824]
# ------------------------------------------------------------------------------
# 1. Load libraries & Set workspace
# ------------------------------------------------------------------------------
library(biomaRt)
library(latex2exp)
library(org.Mm.eg.db)
library(dplyr)
library(GO.db)
library(clusterProfiler)
library(ggplot2)
library(ggrepel)

setwd("/home/local/hoeyoungkim_000504/coursework/")

# ------------------------------------------------------------------------------
# 2. Match annotation to raw counts
# ------------------------------------------------------------------------------
## Load raw counts data generated with FeatureCounts (Multimapping primary)
raw.counts <- read.csv("./SNU-BI1.Data/raw.counts.mm.txt", sep="\t", header=TRUE)

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

## Remove version info from gene IDs
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
rna.counts1 <- meta.counts$RNA.control
rna.counts2 <- meta.counts$RNA.siLin28a
rna.counts3 <- meta.counts$RNA.siLuc

## Filter 2. <80 raw footprint tags in siLuc library
ribo.counts <- meta.counts$RPF.siLuc

filt.meta.counts <- meta.counts[rna.counts1 >= 30 &
                                rna.counts2 >= 30 &
                                rna.counts3 >=30 & 
                                ribo.counts >= 80, ]
nrow(filt.meta.counts)  # 55158 -> 9004

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

norm.density.silin28a <- density.silin28a / 
  (sum(filt.meta.counts$RPF.siLin28a) + sum(filt.meta.counts$RNA.siLin28a))
norm.density.siluc <- density.siluc / 
  (sum(filt.meta.counts$RPF.siLuc) + sum(filt.meta.counts$RNA.siLuc))
rden.change <- log2(norm.density.silin28a / norm.density.siluc)

## Merge data to counts
filt.meta.counts$CLIP.Enrichment <- clip.enrichment
filt.meta.counts$Rden.Change <- rden.change
clip.rden2 <- data.frame(clip.enrichment, rden.change)

## Filter any +/- Inf values
filt.meta.counts <- filt.meta.counts[!is.infinite(filt.meta.counts$CLIP.Enrichment), ]

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
              axis.text.x=element_text(size=11, color="black"),
              axis.text.y=element_text(size=11, color="black"),
              axis.ticks=element_line(size=0.25),
              axis.line=element_line(size=0.25),
              axis.line.x.top=element_blank(),      # Remove top axis line
              axis.line.y.right=element_blank())
ggsave(scatter, filename="./SNU-BI1/FreeExpansionPlots/fig4d.ver3.png", 
       width=4, height=4, units='in', dpi=600)

# ------------------------------------------------------------------------------
# 5. Map GO terms
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
mapping.table <- go.matched[, c(3, 5)]   
length(unique(mapping.table$DB_Object_Symbol)) # Total 7760 / 9004 genes are mapped

## Map GO term to GO description
go.terms <- mapping.table$GO_ID
all.terms <- Term(GOTERM)
matched.descriptions <- all.terms[go.terms]
mapping.table$GO_Description <- matched.descriptions

## Remove any redundant rows
mapping.table <- mapping.table %>% distinct()

# ------------------------------------------------------------------------------
# 6. GO term-to-Genelist Conversion
# ------------------------------------------------------------------------------
go.to.genes <- mapping.table %>% 
  group_by(GO_ID) %>%
  summarise(GeneList = list(unique(DB_Object_Symbol)),
            GO_Description = unique(GO_Description))  
nrow(go.to.genes)   # 14128 unique GO terms  

# ------------------------------------------------------------------------------
# 7. Mann-Whitney U test for each GO term
# ------------------------------------------------------------------------------
## Define a function for Mann-Whitney U test
mw_rden <- function(go_to_genes, rden_df){
  
  ## Initialize a list to store p-values and number of genes
  results <- list()
  num.genes.list <- list()

  for (i in 1:nrow(go_to_genes)) {
    ## Define current GO term and Gene set
    go.term <- go_to_genes$GO_ID[i]
    gene.set <- go_to_genes$GeneList[[i]]
    num.genes <- length(gene.set)
    
    ## Create a logical vector
    rden_df$In_GO <- rden_df$MGIsymbol %in% gene.set
    
    ## Perform the Mann-Whitney U test
    mw.test <- wilcox.test(Rden.Change ~ In_GO, data=rden_df)
    
    # Store the p-value and number of genes
    results[[go.term]] <- mw.test$p.value
    num.genes.list[[go.term]] <- num.genes
  }
  ## # Convert the list of p-values to a vector
  p.values <- unlist(results)
  num.genes.list <- unlist(num.genes.list)
  
  ## Adjust the p-values using the Benjamini-Hochberg method
  p.adj <- p.adjust(p.values, method="BH")
  
  ## Combine results into a dataframe
  results.df <- data.frame(
    GO_ID=names(results),
    GO_Description=go_to_genes$GO_Description,
    Num_Genes=num.genes.list,
    pvalue=p.values,
    p.adj=p.adj
  )
  ## Return results dataframe
  return(results.df)
}

## Define a Rden df (MGIsymbol : Rden.Change)
rden.df <- filt.meta.counts[, c(3, 13)]

## Run
mw.results <- mw_rden(go_to_genes=go.to.genes,
                      rden_df=rden.df)

## Sort according to P.adj
mw.results.sorted <- mw.results[order(mw.results$p.adj), ]
head(mw.results.sorted, 30)

## Save Mann-Whitney U test results
write.table(mw.results.sorted, "./SNU-BI1.Data/mann.whitney.results.csv", sep=",")

# ------------------------------------------------------------------------------
# 8. Filter Mann-Whitney U test data & Remove subset terms
# ------------------------------------------------------------------------------
## Select GO terms with FDR < 0.05
sig.mw.results <- mw.results.sorted[which(abs(mw.results.sorted$p.adj)<0.05),]
nrow(sig.mw.results)   # 446 unique GO terms

## Retrive all Gene sets for each GO terms
sig.go.to.genes <-  go.to.genes[go.to.genes$GO_ID %in% sig.mw.results$GO_ID, ]

## For each GeneList, check if subset
subset <- vector("numeric", nrow(sig.go.to.genes))
subset.of <- vector("numeric", nrow(sig.go.to.genes))

for (i in 1:nrow(sig.go.to.genes)) {
  gene.set <- unlist(sig.go.to.genes$GeneList[i])
  is_subset <- FALSE
  for (j in setdiff(1:nrow(sig.go.to.genes), i)) {
    other.gene.set <- unlist(sig.go.to.genes$GeneList[j])
    if (all(gene.set %in% other.gene.set)) {
      is_subset <- TRUE
      subset.of[i] <- j
      break
    }
  }
  subset[i] <- ifelse(is_subset, 1, 0)
}

sig.go.to.genes$IsSubset <- subset
sig.go.to.genes$SubsetOf <- subset.of

## Remove all subset terms
sig.non.subset <- sig.go.to.genes[sig.go.to.genes$IsSubset==0, ]
nrow(sig.non.subset)  # 446 -> 383 

## Filter M-W test results
sig.mw.non.subset <- sig.mw.results[sig.mw.results$GO_ID %in% sig.non.subset$GO_ID, ]

# ------------------------------------------------------------------------------
# 9. Calculate average CLIP enrichment and Rden Change
# ------------------------------------------------------------------------------
## Average CLIP enrichment for each GO term
clip.enrich.list <- list()
rden.change.list <- list()

for (i in 1:length(sig.mw.non.subset$GO_ID)){
  ## Define current GO term
  go.term <- sig.mw.non.subset$GO_ID[i]
  ## Bring gene sets for current GO term
  gene.set <- unlist(go.to.genes[go.to.genes$GO_ID == go.term, ]$GeneList)
  
  ## Build matched counts-df
  matched <- filt.meta.counts[filt.meta.counts$MGIsymbol %in% gene.set, ]
  
  ## Calculate average values
  avg.clip.enrichment <- sum(matched$CLIP.Enrichment) / nrow(matched)
  avg.rden.change <- sum(matched$Rden.Change) / nrow(matched)
  
  clip.enrich.list[[i]] <- avg.clip.enrichment
  rden.change.list[[i]] <- avg.rden.change
}
clip.enrich.list <- unlist(clip.enrich.list)
rden.change.list <- unlist(rden.change.list)

## Add values to MW dataframe
sig.mw.non.subset$Avg_CLIP <- clip.enrich.list
sig.mw.non.subset$Avg_Rden <- rden.change.list

# ------------------------------------------------------------------------------
# 10. Plot Bubble Plot
# ------------------------------------------------------------------------------
# Reorder by p.adj so smaller p-values are plotted on top
sig.mw.non.subset <- sig.mw.non.subset[order(-sig.mw.non.subset$p.adj), ]

# Define breaks for the color scale in log10 space
log_breaks <- 10^seq(0, -50, by=-5)

# Select the target pathways for labeling
target.GO <- sig.mw.non.subset[c("GO:0005737", "GO:0005634", "GO:0006334",
                                 "GO:0005794", "GO:0005509", "GO:0009986",
                                 "GO:0005788", "GO:0005539", "GO:0005789",
                                 "GO:0005576", "GO:0005783", "GO:0005739"), ]
# Label for target GO
target.label <- list()
for (i in 1:length(target.GO$GO_ID)){
  label <- paste0(target.GO$GO_Description[i], "\n",
                  target.GO$GO_ID[i],
                 " (", target.GO$Num_Genes[i], ")", sep="")
  target.label[i] <- label
}
target.GO$Label <- target.label

## Bubble Plot
bubble <- ggplot(sig.mw.non.subset, aes(x=Avg_CLIP, y=Avg_Rden, size=Num_Genes, colour=p.adj)) +
            geom_point(alpha=0.8) +
            geom_label_repel(data=target.GO, aes(label=Label), size=2.5, segment.size=0.2, 
                            box.padding=1.5, point.padding=1, max.overlaps=Inf,
                            color="#003da8ff", alpha=0.8) +
            coord_cartesian(xlim=c(-2, 2), ylim=c(-2, 1), expand=FALSE) +
            scale_size_area(max_size=15, guide=FALSE) +
            scale_color_gradientn(colors=c("#003da8ff", "#e5eeffff"),
                                  trans="log10",
                                  limits=c(1e-50, 1),
                                  breaks=log_breaks,
                                  labels=scales::trans_format("log10", scales::math_format(10^.x)),
                                  oob=scales::squish) +
            guides(colour=guide_colourbar(reverse=TRUE, barheight=unit(3.5, "inch"),
                                          label.position="right", 
                                          title.position="left", 
                                          title.theme=element_text(angle=90, size=9.5, vjust=0.2))) +
            labs(x=TeX("Enrichment level of LIN28A-bound CLIP tags ($\\log_2$)"),
                 y=TeX("Ribosome density change upon $\\textit{Lin28a}$ knockdown ($\\log_2$)"),
                 color="Term-specific enrichment confidence (False Discovery Rate)") +
            scale_x_continuous(breaks=seq(-2, 2, by=0.5)) +
            scale_y_continuous(breaks=seq(-2, 1, by=0.5)) +
            theme(plot.background=element_rect(fill = "white"),
                  panel.background=element_rect(fill=NA),
                  panel.grid.major=element_line(colour="grey", linetype="dashed", size=0.25),
                  panel.border = element_rect(color = "black", fill = NA, size=0.5),
                  plot.title=element_text(family="Arial", size=9.5),
                  axis.title.x=element_text(family="Arial", size=10),
                  axis.title.y=element_text(family="Arial", size=10),
                  axis.text.x=element_text(size=8, color="black"),
                  axis.text.y=element_text(size=8, color="black"),
                  axis.ticks=element_line(size=0.25),
                  axis.ticks.x.top=element_line(size=0.25),
                  axis.ticks.y.right=element_line(size=0.25),
                  axis.ticks.length=unit(-2, "pt"))

ggsave(bubble, filename="./SNU-BI1/FreeExpansionPlots/bubble.png", 
       width=9, height=4.5, units='in', dpi=300)

