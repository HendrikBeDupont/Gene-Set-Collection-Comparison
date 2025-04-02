# Script:       Final Report.R
# Description:  In this script, we will explore a differential gene 
#               expression dataset comparing disease vs. control
#               samples, perform pathway enrichment analysis and 
#               compare gene set collections.
# Version: 1.0
# Last updated: 2025-02-05
# Author: hdupont


# ##################################################################
# R INSTRUCTIONS
# ##################################################################

# * Lines that start with a # are comments
# * You can run a code line by placing the cursor in the line and clicking 
#   CTRL/Command + Enter

# ##################################################################
# R SETUP
# ##################################################################

# Here we install and load all required packages. 
if (!("BiocManager" %in% installed.packages())) { install.packages("BiocManager", update=FALSE) }
if (!("rstudioapi" %in% installed.packages())) { BiocManager::install("rstudioapi", update=FALSE) }
if (!("org.Hs.eg.db" %in% installed.packages())) { BiocManager::install("org.Hs.eg.db", update=FALSE) }
if (!("dplyr" %in% installed.packages())) { BiocManager::install("dplyr", update=FALSE) }
if (!("EnhancedVolcano" %in% installed.packages())) { BiocManager::install("EnhancedVolcano", update=FALSE) }
if (!("readxl" %in% installed.packages())) { BiocManager::install("readxl", update=FALSE) }
if (!("clusterProfiler" %in% installed.packages())) { BiocManager::install("clusterProfiler", update=FALSE) }
if (!("enrichplot" %in% installed.packages())) { BiocManager::install("enrichplot", update=FALSE) }
if (!("Rgraphviz" %in% installed.packages())) { BiocManager::install("Rgraphviz", update=FALSE) }
if (!("RCy3" %in% installed.packages())) { BiocManager::install("RCy3", update=FALSE) }
if (!("msigdbr" %in% installed.packages())) { BiocManager::install("msigdbr",update=FALSE) }
if (!("RColorBrewer" %in% installed.packages())) { BiocManager::install("RColorBrewer",update=FALSE) }
if (!("readr" %in% installed.packages())) { BiocManager::install("readr",update=FALSE) }
if (!("rWikiPathways" %in% installed.packages())) { BiocManager::install("rWikiPathways",update=FALSE) }
if (!("ReactomePA" %in% installed.packages())) { BiocManager::install("ReactomePA", update=FALSE) }

library(rstudioapi)
library(org.Hs.eg.db)
library(dplyr)
library(EnhancedVolcano)
library(readxl)
library(clusterProfiler)
library(enrichplot)
library(Rgraphviz)
library(RCy3)
library(msigdbr)
library(RColorBrewer)
library(readr)
library(ggplot2)
library(rWikiPathways)
library(ReactomePA)

# We will set the working directory to the location where the current 
# script is located. This way, we can use relative file path locations. 
setwd(getwd())

# We will create an output folder where all figures and files will be stored
out.folder <- "output/"
dir.create(out.folder)

# ##################################################################
# 1.1 Data Import and Preprocessing
# ##################################################################

# Load dataset
data <- read_excel("GSE239914-differential-analysis.xlsx")

# Select specific columns (gene information, logFC, p-value)
data <- data[,c(8,1,6,2,3,10,11)]

# Filter for protein-coding genes and remove missing values
data.pc <- data %>%
  filter(GeneType == "protein-coding") %>%
  filter(!is.na(GeneID))

# Dropping a column 6 from the dataset 
data.pc <- data.pc[,-6]

# ##################################################################
# Decision: Choose log2FC Cutoff, Set confidence Cutoff
# ##################################################################

# Define thresholds for differential expression
log2fc.cutoff <- 1
pvalue.cutoff <- 0.05

# ##################################################################
# 1.2 Filtering for Differentially Expressed Genes (DEGs)
# ##################################################################

# Select differently expressed genes (DEGs) (genes meeting log2FC and p-value criteria)
degs <- data.pc[abs(data.pc$log2FoldChange) > log2fc.cutoff & data.pc$padj < pvalue.cutoff,]

# Save the DEGs to a file
write.table(degs, file=paste0(out.folder,"degs.tsv"), row.names = FALSE, sep="\t", quote = FALSE)

# ##################################################################
# Optional: Run analysis on up- or down-regulated genes
# ##################################################################

# Selects genes where log2FoldChange is greater than the threshold (1 by default)
genes.up <- nrow(degs[degs$log2FoldChange > log2fc.cutoff,])

# Selects genes where log2FoldChange is less than -1
genes.down <- nrow(degs[degs$log2FoldChange < -log2fc.cutoff,])

# Print results
cat("Number of Upregulated Genes: ", genes.up, "\n")
cat("Number of Downregulated Genes: ", genes.down, "\n")

# ##################################################################
# 2.1 Over-representation analysis (ORA)
# ##################################################################

# Fetch WikiPathways gene sets
genesets.wp <- msigdbr(species = "Homo sapiens", subcategory = "CP:WIKIPATHWAYS") %>% dplyr::select(gs_name, entrez_gene)

# Perform pathway enrichment analysis
res.wp <- clusterProfiler::enricher(degs$GeneID, TERM2GENE = genesets.wp, pAdjustMethod = "fdr", pvalueCutoff = 0.05, minGSSize = 5, maxGSSize = 400)
res.wp.df <- as.data.frame(res.wp)

# Save results 
write.table(res.wp.df, file=paste0(out.folder,"WP-Enrichment.txt"), sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# ##################################################################
# 2.2 Run enrichment for each gene set collection
# ##################################################################

# GO Enrichment (Biological Process)
ego <- enrichGO(gene = degs$GeneID,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "BP",
                pAdjustMethod = "fdr",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego.df <- as.data.frame(ego)
write.table(ego.df, file=paste0(out.folder,"GO-Enrichment.txt"), sep="\t", row.names=FALSE, quote=FALSE)

# KEGG Enrichment
ekegg <- enrichKEGG(gene = degs$GeneID,
                    organism = "hsa",
                    pvalueCutoff = 0.05)
ekegg.df <- as.data.frame(ekegg)
write.table(ekegg.df, file=paste0(out.folder,"KEGG-Enrichment.txt"), sep="\t", row.names=FALSE, quote=FALSE)

# MSigDB Hallmark Gene Sets (H)
msig.h <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, entrez_gene)
emsig <- enricher(degs$GeneID,
                  TERM2GENE = msig.h,
                  pAdjustMethod = "fdr",
                  pvalueCutoff = 0.05)
emsig.df <- as.data.frame(emsig)
write.table(emsig.df, file=paste0(out.folder,"MSigDB-Enrichment.txt"), sep="\t", row.names=FALSE, quote=FALSE)

# Reactome Enrichment
ereactome <- enrichPathway(gene = degs$GeneID,
                           organism = "human",
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "fdr",
                           minGSSize = 5,
                           readable = TRUE)
ereactome.df <- as.data.frame(ereactome)
write.table(ereactome.df, file = paste0(out.folder, "Reactome-Enrichment.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

# ##################################################################
# 2.3 Count the number of significantly enriched terms
# ##################################################################

n.go <- nrow(ego.df)
n.kegg <- nrow(ekegg.df)
n.wp <- nrow(res.wp.df)
n.msig <- nrow(emsig.df)
n.reactome <- nrow(ereactome.df)

# Combine into one dataframe
enrichment_counts <- data.frame(
  Database = c("GO (BP)", "KEGG", "WikiPathways", "MSigDB (Hallmark)", "Reactome"),
  SignificantTerms = c(n.go, n.kegg, n.wp, n.msig, n.reactome)
)

# Save the table
write.table(enrichment_counts, file=paste0(out.folder,"Enrichment_Comparison.txt"), sep="\t", row.names=FALSE, quote=FALSE)

# See table
print(enrichment_counts)

# ##################################################################
# 2.4 Visualize the comparison with a barplot
# ##################################################################

# Barplot to visualize the comparison
ggplot(enrichment_counts, aes(x = Database, y = SignificantTerms, fill = Database)) +
  geom_bar(stat = "identity") +
  labs(title = "Comparison of Significantly Enriched Terms",
       y = "Number of Significant Pathways",
       x = "") +
  theme_minimal() +
  theme(legend.position = "none")

# Save plot
ggsave(filename = paste0(out.folder, "Enrichment_Comparison_Barplot.png"), width = 8, height = 6)

# ##################################################################
# 3.1 Pathway Relevance Scoring
# ##################################################################

# Define vEDS-relevant biological keywords
keywords <- c(
  # Core structural components
  "collagen", "extracellular matrix", "fibril", "fibrillar", "glycoprotein", "elastin", 
  "basement membrane", "matrisome",
  
  # Cellular structures and adhesion
  "focal adhesion", "cell adhesion", "integrin", "junction", "cytoskeleton",
  
  # Tissue types and general biology
  "connective tissue", "smooth muscle", "vascular", "blood vessel", "angiogenesis", 
  "endothelial", "pericyte", "myofibroblast",
  
  # Pathophysiological processes
  "wound healing", "inflammation", "fibrosis", "remodeling", "matrix degradation", 
  "proteoglycan", "MMP", "TIMP", "TGF", "ECM disassembly",
  
  # Genetic/rare disease links
  "Ehlers", "Marfan", "cutis laxa", "arterial dissection", "aneurysm", 
  "vascular development", "vascular integrity"
)

# ##################################################################
# 3.2 Search for relevant pathways in each result
# ##################################################################

# Prepare your enrichment result data frames
# 2.1 + A.1 

# Function to count relevant terms based on keyword list
count_relevant_terms <- function(df, term_column = "Description", keywords) {
  keyword_matches <- sapply(keywords, function(kw) {
    grepl(kw, tolower(df[[term_column]]), ignore.case = TRUE)
  })
  relevance_score <- rowSums(keyword_matches) > 0
  return(sum(relevance_score))
}

# Count relevant pathways
relevance_scores <- data.frame(
  Database = c("GO", "KEGG", "WikiPathways", "MSigDB", "Reactome"),
  RelevantTerms = c(
    count_relevant_terms(ego.df, "Description", keywords),
    count_relevant_terms(ekegg.df, "Description", keywords),
    count_relevant_terms(res.wp.df, "Description", keywords),
    count_relevant_terms(emsig.df, "Description", keywords),
    count_relevant_terms(ereactome.df, "Description", keywords)
  )
)

# Save the table
write.table(relevance_scores, file = paste0(out.folder, "Pathway_Relevance_Score.txt"), sep = "\t", row.names = FALSE)

# ##################################################################
# 3.3 Visualize the comparison
# ##################################################################

# Barplot to compare pathway relevance scores
ggplot(relevance_scores, aes(x = Database, y = RelevantTerms, fill = Database)) +
  geom_bar(stat = "identity") +
  labs(title = "Pathway Relevance Scoring (vEDS-related Pathways)",
       y = "Number of vEDS-Relevant Pathways",
       x = "") +
  theme_minimal() +
  theme(legend.position = "none")

# Save plot
ggsave(filename = paste0(out.folder, "Pathway_Relevance_Scoring.png"), width = 8, height = 6)

# ##################################################################
# 3.4 Comparison Table
# ##################################################################

# Function to extract most significant p-adjusted value among relevant terms
get_most_significant_padj <- function(df, keywords) {
  matches <- apply(df, 1, function(row) {
    any(sapply(keywords, function(kw) grepl(kw, tolower(row["Description"]), ignore.case = TRUE)))
  })
  matched <- df[matches, ]
  if (nrow(matched) == 0) return(NA)
  return(formatC(min(matched$p.adjust, na.rm = TRUE), format = "e", digits = 2))
}

# Function to get the top relevant term name (lowest padj)
get_top_term <- function(df, keywords) {
  matches <- apply(df, 1, function(row) {
    any(sapply(keywords, function(kw) grepl(kw, tolower(row["Description"]), ignore.case = TRUE)))
  })
  matched <- df[matches, ]
  if (nrow(matched) == 0) return(NA)
  top_row <- matched[which.min(matched$p.adjust), ]
  return(top_row$Description)
}

# Create the complete comparison table
final_table <- data.frame(
  Collection = c("GO", "KEGG", "WikiPathways", "MSigDB", "Reactome"),
  ECM_Related_Terms_Found = relevance_scores$RelevantTerms,
  Most_Significant_padj = c(
    get_most_significant_padj(ego.df, keywords),
    get_most_significant_padj(ekegg.df, keywords),
    get_most_significant_padj(res.wp.df, keywords),
    get_most_significant_padj(emsig.df, keywords),
    get_most_significant_padj(ereactome.df, keywords)
  ),
  Top_Term = c(
    get_top_term(ego.df, keywords),
    get_top_term(ekegg.df, keywords),
    get_top_term(res.wp.df, keywords),
    get_top_term(emsig.df, keywords),
    get_top_term(ereactome.df, keywords)
  )
)

# Save the final table
write.table(final_table, file = paste0(out.folder, "ECM_Relevance_Table.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

# View it for copy/paste
View(final_table)

# ##################################################################
# 4.1 Functional Clustering for GO (with simplify)
# ##################################################################

# Simplify GO enrichment result to reduce redundancy
ego.simplified <- simplify(
  ego,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang"
)

# Compute similarity
ego.sim <- pairwise_termsim(ego.simplified)

# Tree plot
treeplot(ego.sim, showCategory = 20, label_format = 50)
ggsave(paste0(out.folder, "GO_Treeplot.png"), width = 10, height = 10)

# Enrichment map
emapplot(ego.sim, showCategory = 20, layout = "kk", cex_category = 0.8)
ggsave(paste0(out.folder, "GO_EnrichmentMap.png"), width = 10, height = 8)

# ##################################################################
# 4.2 KEGG Clustering
# ##################################################################

ekegg.sim <- pairwise_termsim(ekegg)

treeplot(ekegg.sim, showCategory = 20, label_format = 50)
ggsave(paste0(out.folder, "KEGG_Treeplot.png"), width = 10, height = 10)

emapplot(ekegg.sim, showCategory = 20, layout = "kk", cex_category = 0.8)
ggsave(paste0(out.folder, "KEGG_EnrichmentMap.png"), width = 10, height = 8)

# ##################################################################
# 4.3 WikiPathways Clustering
# ##################################################################

res.wp.sim <- pairwise_termsim(res.wp)

treeplot(res.wp.sim, showCategory = 20, label_format = 50)
ggsave(paste0(out.folder, "WikiPathways_Treeplot.png"), width = 10, height = 10)

emapplot(res.wp.sim, showCategory = 20, layout = "kk", cex_category = 0.8)
ggsave(paste0(out.folder, "WikiPathways_EnrichmentMap.png"), width = 10, height = 8)

# ##################################################################
# 4.4 MSigDB Clustering
# ##################################################################

emsig.sim <- pairwise_termsim(emsig)

treeplot(emsig.sim, showCategory = 20, label_format = 50)
ggsave(paste0(out.folder, "MSigDB_Treeplot.png"), width = 10, height = 10)

emapplot(emsig.sim, showCategory = 20, layout = "kk", cex_category = 0.8)
ggsave(paste0(out.folder, "MSigDB_EnrichmentMap.png"), width = 10, height = 8)

# ##################################################################
# 4.5 Reactome Clustering
# ##################################################################

ereactome.sim <- pairwise_termsim(ereactome)

treeplot(ereactome.sim, showCategory = 20, label_format = 50)
ggsave(paste0(out.folder, "Reactome_Treeplot.png"), width = 10, height = 10)

emapplot(ereactome.sim, showCategory = 20, layout = "kk", cex_category = 0.8)
ggsave(paste0(out.folder, "Reactome_EnrichmentMap.png"), width = 10, height = 8)


# ================================================
# FUNCTION: functional_clustering_plots
# ================================================
functional_clustering_plots <- function(enrich_result, database_name, out.folder, show_n = 20, simplify_go = FALSE) {
  # Optional simplification for GO
  if (simplify_go) {
    enrich_result <- simplify(
      enrich_result,
      cutoff = 0.7,
      by = "p.adjust",
      select_fun = min,
      measure = "Wang",
      ontology = "BP"
    )
  }
  
  # Compute similarity
  enrich_sim <- enrichplot::pairwise_termsim(enrich_result)
  
  # Tree plot
  treeplot_path <- paste0(out.folder, database_name, "_Treeplot.png")
  png(treeplot_path, width = 1000, height = 1000)
  print(treeplot(enrich_sim, showCategory = show_n, label_format = 50))
  dev.off()
  
  # Enrichment map
  emapplot_path <- paste0(out.folder, database_name, "_EnrichmentMap.png")
  png(emapplot_path, width = 1000, height = 800)
  print(emapplot(enrich_sim, showCategory = show_n, layout = "kk", cex_category = 0.8))
  dev.off()
  
  message(paste("Plots saved for", database_name))
}

# example usage for all databses
functional_clustering_plots(ego, "GO", out.folder, simplify_go = TRUE)
functional_clustering_plots(ekegg, "KEGG", out.folder)
functional_clustering_plots(res.wp, "WikiPathways", out.folder)
functional_clustering_plots(emsig, "MSigDB", out.folder)
functional_clustering_plots(ereactome, "Reactome", out.folder)
