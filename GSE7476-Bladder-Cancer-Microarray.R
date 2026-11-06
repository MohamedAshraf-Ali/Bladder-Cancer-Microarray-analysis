#Loading required packages
library(affy)
library(affyPLM)
library(Biobase)
library(readxl)
library(tidyverse)
library(matrixStats)
library(ggfortify)
library(limma)
library(AnnotationDbi)
library(hgu133plus2.db)
library(ComplexHeatmap)
library(EnhancedVolcano)
library(plotly)
library(ComplexHeatmap)
library(circlize)
library(limma)
library(ggVennDiagram)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(msigdbr)


#Define input directories 
cel_dir <- "GSE7476"
pheno_file <- "phenotable.xlsx"

#check input files 
stopifnot(dir.exists(cel_dir))
stopifnot(file.exists(pheno_file))

#Load CEL files 
affy_data <- ReadAffy( celfile.path = cel_dir)

#Load meta data
pheno_data <- as.data.frame(read_excel(pheno_file))
pheno_data$Sample.Type <- unlist(sapply(strsplit(pheno_data$Sample.Type, "_"), c) [2, ])

#Inspect metadata 
head(pheno_data)
str(pheno_data)
colnames(pheno_data)

#Inspect affy_data
class(affy_data)
dim(affy_data)
sampleNames(affy_data)
featureNames(affy_data)
annotation(affy_data)
View(exprs(affy_data))

#Check the number of CEL files
length(sampleNames(affy_data))

#Check duplicated sample names 
anyDuplicated(sampleNames(affy_data))

#Check sample group distribution 
table(pheno_data$Sample.Type)

#Check samples missing from both: phenotable and CEL data
# Standardize CEL sample names
cel_sample_ids <- sub("\\.CEL\\.gz$", "", sampleNames(affy_data), ignore.case = TRUE)
setdiff( cel_sample_ids, pheno_data$Sample.Id)
setdiff(pheno_data$Sample.Id, cel_sample_ids)

#Ensure phenotype data matches CEL files 
stopifnot(identical(cel_sample_ids, pheno_data$Sample.Id))

################## Raw Microarray Quality Control ####################
exp.raw <- Biobase::exprs(affy_data)

#Raw Density plot  
raw_plot_data <- as.data.frame(exp.raw) %>%
  rownames_to_column("Probe_IDs") %>%
  pivot_longer(
    cols = -Probe_IDs, names_to = "Sample IDs", values_to = "Intensity"
  )

ggplot(
  raw_plot_data, aes(x = Intensity, colour = `Sample IDs`)
) + 
  geom_density( linewidth = 0.7) + 
  labs(
    title = "Raw Microarray Intensity Distribution",
    x = "Probe Intensity", y = "Density", colour = "Samples"
  ) + theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  filename = "Raw_Density_Plot.png",
  plot = raw_density_plot,
  width = 10,
  height = 7,
  units = "in"
)

#Intensity BoxPlot 
ggplot(raw_plot_data, aes(x = `Sample IDs`, y = Intensity, fill = `Sample IDs`)) +
  geom_boxplot(outlier.size = 0.3, linewidth = 0.3) +
  labs(
    title = "Raw Microarray Intenstiy Distribution across Samples",
    x = "Sample IDs", y = "Probe Intensity"
  ) + theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  )

ggsave(
  filename = "Raw_Boxplot.png",
  plot = raw_boxplot,
  width = 10,
  height = 7
)

#PCA 
pca_raw <- prcomp(t(exp.raw), scale. = F, center = T)
autoplot(pca_raw, data = pheno_data, colour = "Sample.Type", frame = T)

ggsave(
  "Raw_PCA.png",
  plot = raw_pca_plot,
  width = 9,
  height = 7
)

#3D PCA 
colnames(exp.raw) <- cel_sample_ids
exp.raw.log2 <- log2(exp.raw + 1)

#Percentage of variance explained by each PCA 
variance_explained <- (pca_raw$sdev^2 / sum(pca_raw$sdev^2)) * 100
variance_explained[1:3]

#Prepare PCA data for plotting 
pca_data <- as.data.frame(pca_raw$x[, 1:3])
pca_data$Sample.IDs <- rownames(pca_data)

pca_data$Sample.Type <- pheno_data$Sample.Type [ 
  match (pca_data$Sample.IDs, pheno_data$Sample.IDs)]

#Plotting 
plot_ly(
  data = pca_data, x = ~ PC1, y = ~ PC2, z = ~ PC3, color = ~Sample.Type, 
  colors = c("#0072B2", "#D55E00"), mode = "markers+text", text = ~Sample.IDs
) %>%
  layout(
    title = list(text = "3D PCA - Raw Affymetrix Expression Data", y = 0.95)
  )

################# Background Correction and Normalization ###################
threestep_data <- affyPLM::threestep( object = affy_data,
  background.method = "IdealMM",
  normalize.method = "quantile",
  summary.method = "average.log"
)

#Extract normalized matrix 
exp_norm <- exprs(threestep_data)
colnames(exp_norm) <- unlist(sapply(strsplit(colnames(exp_norm), "\\."), c) [1,])

#Density plot for normalized matrix 
norm_plot_data <- as.data.frame(exp_norm) %>%
  rownames_to_column("Probe_IDs") %>%
  pivot_longer(
    cols = -Probe_IDs, names_to = "Sample.IDs", values_to = "Expression"
  )

ggplot(
  norm_plot_data, aes(x = Expression, colour = Sample.IDs)
) + geom_density(linewidth = 0.7) +
  labs(
    title = "Three-Step Normalized Expression Distribution",
    x = "Expression intensity",
    y = "Density",
    colour = "Samples"
  ) + theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  filename = "Normalized_Density_Plot.png",
  plot = norm_density_plot,
  width = 10,
  height = 7,
  units = "in"
)

#Boxplot for normalized matrix 
ggplot(
  norm_plot_data, aes(x = Sample.IDs, y = Expression, fill = Sample.IDs)
) + geom_boxplot(
  linewidth = 0.3, outlier.size = 0.3
) + labs(
  title = "Three-Step Normalized Expression Values",
  x = "Samples",
  y = "Expression intensity",
) + theme(
  plot.title = element_text(hjust = 0.5),
  axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 0.5
  )
)

ggsave(
  "Normalized_Boxplot.png",
  plot = norm_boxplot,
  width = 10,
  height = 7
)

#Sample correlation heatmap
#Calculate Pearson correlation between samples 
sample_cor <- cor( exp_norm, method = "pearson", use = "pairwise.complete.obs")
round(sample_cor, 3)

# Match sample groups to expression matrix columns
sample_annotation <- pheno_data [match(colnames(exp_norm), pheno_data$Sample.Id), ]
rownames(sample_annotation) <- sample_annotation$Sample.Id

#Create the correlation heatmap 
#Define group colors 
group_colors <- c ("TaGL"   = "#0072B2",
                   "T1GH"   = "#D55E00",
                   "T2-4"   = "#009E73",
                   "control" = "#CC79A7")

#Create column annotation 
ha <- HeatmapAnnotation(Sample.Type = sample_annotation$Sample.Type, 
                        col = list(Sample.Type = group_colors))
#Correlation Heatmap 
Heatmap(sample_cor, name = "Pearson\nCorrelation",
        col = colorRamp2(c(0.5, 0.75, 1), c("#2166AC", "#F7F7F7", "#B2182B")),
        top_annotation = ha, cluster_rows = T, cluster_columns = T,
        column_title = "Sample-to-Sample Correlation", row_title = "Columns",
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))

png(
  "Sample_Correlation_Heatmap.png",
  width = 2700,
  height = 2400,
  res = 300
)

draw(correlation_heatmap)
dev.off()

################### Probe Annotation #####################
#Extract probe IDs 
probe_IDs <- rownames(exp_norm)

#Retrieve annotations 
gene_annotation <- AnnotationDbi::select(hgu133plus2.db, keys = probe_IDs,
                                         columns = c("SYMBOL", "GENENAME", "ENTREZID"),
                                         keytype = "PROBEID")

#Check annotated coverage 
length(probe_IDs)
sum(is.na(gene_annotation$SYMBOL))

#Remove probes with missing gene annotations 
gene_annotation_clean <- gene_annotation %>%
  filter(!is.na(SYMBOL), SYMBOL != "")

#Check duplicated probe IDs 
sum(duplicated(gene_annotation_clean$PROBEID))

# Inspect duplicated records
gene_annotation_clean %>%
  dplyr::count(PROBEID) %>%
  dplyr::filter(n > 1) %>%
  dplyr::arrange(desc(n))

# Number of probes in expression matrix
nrow(exp_norm)

# Number of unique annotated probes
length(unique(gene_annotation_clean$PROBEID))

# Number of probes without valid gene symbols
length(setdiff(
  rownames(exp_norm),
  unique(gene_annotation_clean$PROBEID)
))

# Keep one annotation record per probe
gene_annotation_unique <- gene_annotation_clean %>%
  dplyr::distinct(PROBEID, .keep_all = TRUE)

# Check the result
nrow(gene_annotation_unique)
length(unique(gene_annotation_unique$PROBEID))

# Confirm that probe IDs are unique
anyDuplicated(gene_annotation_unique$PROBEID)

# Check which probes are present in both objects
sum(rownames(exp_norm) %in% gene_annotation_unique$PROBEID)

# Check probes missing from the annotation table
sum(!rownames(exp_norm) %in% gene_annotation_unique$PROBEID)

# Check annotation probes absent from the expression matrix
sum(!gene_annotation_unique$PROBEID %in% rownames(exp_norm))

gene_annotation_matched <- gene_annotation_unique [match(rownames(exp_norm), 
                                                         gene_annotation_unique$PROBEID), ]

# Verify alignment
stopifnot(
  identical(
    rownames(exp_norm),
    gene_annotation_matched$PROBEID
  )
)

# Check the number of rows
nrow(exp_norm)
nrow(gene_annotation_matched)

# Restore all original probe IDs, including unannotated probes
gene_annotation_matched$PROBEID <- rownames(exp_norm)

# Verify alignment
stopifnot(
  identical(rownames(exp_norm), gene_annotation_matched$PROBEID)
)

#Prepare the expression matrix for differential expression analysis
#Define the sample group as a factor 
pheno_data$Sample.Type <- factor(pheno_data$Sample.Type, 
                                 levels = c("control", "TaGL", "T1GH", "T2-4"))

#Build the design matrix 
design <- model.matrix( ~0 + Sample.Type, data = pheno_data)
colnames(design) <- levels(pheno_data$Sample.Type)

#Fit the linear model: limma package 
fit <- limma::lmFit(object = exp_norm, design = design)

# Rename design matrix columns
colnames(design) <- c(
  "control",
  "TaGL",
  "T1GH",
  "T2_4"
)

 #Define differential expression comparison: The three tumour groups vs the control group
contrasts <- limma::makeContrasts(
  TaGL_vs_control = TaGL - control,
  T1GH_vs_control = T1GH - control,
  T2_4_vs_control = T2_4 - control,
  levels = design
)

#Apply contrasts to the fitted limma model
fit_contrasts <- limma::contrasts.fit(fit = fit, contrasts = contrasts)

#Apply empirical Bayes moderation 
fit_ebayes <- limma::eBayes(fit_contrasts)

# Inspect the fitted model
summary(fit_ebayes)

# Check the number of significantly differentially expressed probes
limma::decideTests(fit_ebayes)

#################### Extract the differential expression results #######################
#TaGL vs control 
deg_TaGL <- limma::topTable(fit_ebayes, coef = "TaGL_vs_control", 
                            number = Inf, adjust.method = "BH", sort.by = "P")

#T1GH vs control 
deg_T1GH <- limma::topTable(fit_ebayes, coef = "T1GH_vs_control",
                            number = Inf, adjust.method = "BH", sort.by = "P")

#T2_4 vs control 
deg_T2_4 <- limma::topTable(fit_ebayes, coef = "T2_4_vs_control",
                            number = Inf, adjust.method = "BH", sort.by = "P")

#Count significant probes using an adjusted p-value threshold of 0.05
sum(deg_TaGL$adj.P.Val < 0.05, na.rm = TRUE)
sum(deg_T1GH$adj.P.Val < 0.05, na.rm = TRUE)
sum(deg_T2_4$adj.P.Val < 0.05, na.rm = TRUE)

#Add gene annotations to each DEG table
deg_TaGL$PROBEID <- rownames(deg_TaGL)
deg_T1GH$PROBEID <- rownames(deg_T1GH)
deg_T2_4$PROBEID <- rownames(deg_T2_4)

deg_TaGL_annotated <- dplyr::left_join(deg_TaGL, gene_annotation_matched, by = "PROBEID")
deg_T1GH_annotated <- dplyr::left_join(deg_T1GH, gene_annotation_matched, by = "PROBEID")
deg_T2_4_annotated <- dplyr::left_join(deg_T2_4, gene_annotation_matched, by = "PROBEID")

################# identify the upregulated and downregulated probes ################
# Define DEG classification function
classify_DEGs <- function(data) {
  dplyr::mutate(
    data,
    Regulation = dplyr::case_when(
      adj.P.Val < 0.05 & logFC > 0  ~ "Upregulated",
      adj.P.Val < 0.05 & logFC < 0 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )
}

# Classify each comparison
deg_TaGL_annotated <- classify_DEGs(deg_TaGL_annotated)
deg_T1GH_annotated <- classify_DEGs(deg_T1GH_annotated)
deg_T2_4_annotated <- classify_DEGs(deg_T2_4_annotated)

#Count the results 
# TaGL vs control
table(deg_TaGL_annotated$Regulation)

# T1GH vs control
table(deg_T1GH_annotated$Regulation)

# T2-4 vs control
table(deg_T2_4_annotated$Regulation)

#Summarize DEG counts 
count_DEGs_by_threshold <- function(data) {
  significant <- data$adj.P.Val < 0.05
  
  c(
    Significant = sum(significant, na.rm = TRUE),
    Upregulated = sum(
      significant & data$logFC > 0,
      na.rm = TRUE
    ),
    Downregulated = sum(
      significant & data$logFC < 0,
      na.rm = TRUE
    ),
    Up_FC_gt_1 = sum(
      significant & data$logFC > 1,
      na.rm = TRUE
    ),
    Down_FC_lt_minus_1 = sum(
      significant & data$logFC < -1,
      na.rm = TRUE
    )
  )
}

deg_summary <- rbind(
  TaGL = count_DEGs_by_threshold(deg_TaGL),
  T1GH = count_DEGs_by_threshold(deg_T1GH),
  T2_4 = count_DEGs_by_threshold(deg_T2_4)
)

#################### Export Annotation and DEG Results ####################

# 1. Probe annotation
write.csv(
  gene_annotation_unique,
  "Gene_Annotation_Unique.csv",
  row.names = FALSE
)

write.csv(
  gene_annotation_matched,
  "Gene_Annotation_Matched.csv",
  row.names = FALSE
)

# 2. Complete annotated differential expression results
write.csv(
  deg_TaGL_annotated,
  "DEG_TaGL_vs_Control.csv",
  row.names = FALSE
)

write.csv(
  deg_T1GH_annotated,
  "DEG_T1GH_vs_Control.csv",
  row.names = FALSE
)

write.csv(
  deg_T2_4_annotated,
  "DEG_T2_4_vs_Control.csv",
  row.names = FALSE
)

# 3. DEG summary
deg_summary_df <- as.data.frame(deg_summary) %>%
  tibble::rownames_to_column("Comparison")

write.csv(
  deg_summary_df,
  "DEG_Summary.csv",
  row.names = FALSE
)

# Verify exported files
list.files(
  pattern = "^(Gene_Annotation|DEG_).*\\.csv$"
)

################ Volcano Plot #################
#Volcano plot for TaGL
volcano_TaGL <- EnhancedVolcano(
  deg_TaGL_annotated,
  lab = deg_TaGL_annotated$SYMBOL,
  x = "logFC",
  y = "adj.P.Val",
  title = "TaGL vs Control",
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 1.5,
  labSize = 3,
  drawConnectors = TRUE,
  max.overlaps = 20
)

ggsave("Volcano_TaGL_vs_Control.png", volcano_TaGL, width = 9, height = 8)

#Volcano plot for T1GH
volcano_T1GH <- EnhancedVolcano(
  deg_T1GH_annotated,
  lab = deg_T1GH_annotated$SYMBOL,
  x = "logFC",
  y = "adj.P.Val",
  title = "T1GH vs Control",
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 1.5,
  labSize = 3,
  drawConnectors = TRUE,
  max.overlaps = 20
)

ggsave("Volcano_T1GH_vs_Control.png", volcano_T1GH, width = 9, height = 8)

#Volcano plot for T2_4
volcano_T2_4 <- EnhancedVolcano(
  deg_T2_4_annotated,
  lab = deg_T2_4_annotated$SYMBOL,
  x = "logFC",
  y = "adj.P.Val",
  title = "T2_4 vs Control",
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 1.5,
  labSize = 3,
  drawConnectors = TRUE,
  max.overlaps = 20
)

ggsave("Volcano_T2_4_vs_Control.png", volcano_T2_4, width = 9, height = 8)


#Identify significant probes from each comparison 
sig_TaGL <- rownames(deg_TaGL) [deg_TaGL$adj.P.Val < 0.05]
sig_T1GH <- rownames(deg_T1GH) [deg_T1GH$adj.P.Val < 0.05]
sig_T2_4 <- rownames(deg_T2_4) [deg_T2_4$adj.P.Val < 0.05]

#Identify shared and unique probes 
#Shared probes across all comparisons 
common_DEGs <- intersect(intersect(sig_TaGL, sig_T1GH), sig_T2_4)

#Probes unique to each comparison
unique_TaGL <- setdiff(sig_TaGL, union(sig_T1GH, sig_T2_4))
unique_T1GH <- setdiff(sig_T1GH, union(sig_TaGL, sig_T2_4))
unique_T2_4 <- setdiff(sig_T2_4, union(sig_T1GH, sig_TaGL))

#Summary of the overlap 
overlap_summary <- data.frame(
  Category = c(
    "TaGL",
    "T1GH",
    "T2-4",
    "Shared across all three",
    "Unique to TaGL",
    "Unique to T1GH",
    "Unique to T2-4"
  ),
  Number_of_Probes = c(
    length(sig_TaGL),
    length(sig_T1GH),
    length(sig_T2_4),
    length(common_DEGs),
    length(unique_TaGL),
    length(unique_T1GH),
    length(unique_T2_4)
  )
)

#Visualize the overlap: Venn Diagram
deg_sets <- list(TaGL = sig_TaGL, T1GH = sig_T1GH, T2_4 = sig_T2_4)
venn_plot <- ggVennDiagram(deg_sets, label_alpha = 0.7) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_void()

ggsave("Significant_Probe_Overlap_Venn.png",
       plot = venn_plot, width = 8, height = 7)

#Do shared DEGs change in the same direction across all three comparisons?
shared_DEG_results <- data.frame(Probe_IDs = common_DEGs,
                                 TaGL_logFC = deg_TaGL [common_DEGs, "logFC"],
                                 T1GH_logFC = deg_T1GH [common_DEGs, "logFC"],
                                 T2_4_logFC = deg_T2_4 [common_DEGs, "logFC"])
shared_DEG_results$Direction <- apply(
  shared_DEG_results[, 2:4],
  1,
  function(x) {
    if (all(x > 0)) {
      "Consistently upregulated"
    } else if (all(x < 0)) {
      "Consistently downregulated"
    } else {
      "Opposite directions"
    }
  }
)

# Save overlap summary
write.csv(
  overlap_summary,
  "Significant_Probe_Overlap_Summary.csv",
  row.names = FALSE
)

# Save shared DEG results
write.csv(
  shared_DEG_results,
  "Shared_DEG_Probe_Results.csv",
  row.names = FALSE
)

# Save individual and shared probe lists
write.csv(data.frame(Probe_ID = sig_TaGL),
          "Significant_Probes_TaGL.csv", row.names = FALSE)

write.csv(data.frame(Probe_ID = sig_T1GH),
          "Significant_Probes_T1GH.csv", row.names = FALSE)

write.csv(data.frame(Probe_ID = sig_T2_4),
          "Significant_Probes_T2_4.csv", row.names = FALSE)

write.csv(data.frame(Probe_ID = common_DEGs),
          "Shared_Significant_Probes.csv", row.names = FALSE)

############ gene-level annotation and deduplication before functional enrichment ############
annotation_summary <- gene_annotation_unique %>%
  summarise(
    Total_probes = n(),
    Unique_genes = n_distinct(SYMBOL),
    Missing_symbols = sum(is.na(SYMBOL) | SYMBOL == "")
  )

#Check how many unique genes are represented in each significant DEG set
count_genes <- function(probes, annotation_table) {
  annotation_table %>%
    dplyr::filter(
      PROBEID %in% probes,
      !is.na(SYMBOL),
      SYMBOL != ""
    ) %>%
    dplyr::summarise(
      Significant_probes = dplyr::n_distinct(PROBEID),
      Unique_genes = dplyr::n_distinct(SYMBOL)
    )
}

count_genes(sig_TaGL, gene_annotation_unique)
count_genes(sig_T1GH, gene_annotation_unique)
count_genes(sig_T2_4, gene_annotation_unique)

#Create gene-level DEG tables
create_gene_table <- function(deg_table, annotation_table) {
  
  deg_table$PROBEID <- rownames(deg_table)
  
  gene_table <- deg_table %>%
    dplyr::inner_join(
      annotation_table %>%
        dplyr::filter(!is.na(SYMBOL), SYMBOL != "") %>%
        dplyr::select(PROBEID, SYMBOL, GENENAME, ENTREZID),
      by = "PROBEID"
    ) %>%
    dplyr::filter(!is.na(adj.P.Val)) %>%
    dplyr::arrange(adj.P.Val) %>%
    dplyr::distinct(SYMBOL, .keep_all = TRUE)
  
  gene_table
}

gene_TaGL <- create_gene_table(deg_TaGL, gene_annotation_unique)
gene_T1GH <- create_gene_table(deg_T1GH, gene_annotation_unique)
gene_T2_4 <- create_gene_table(deg_T2_4, gene_annotation_unique)

#Significant gene-level tables using FDR < 0.05
sig_gene_TaGL <- gene_TaGL %>%
  dplyr::filter(adj.P.Val < 0.05)

sig_gene_T1GH <- gene_T1GH %>%
  dplyr::filter(adj.P.Val < 0.05)

sig_gene_T2_4 <- gene_T2_4 %>%
  dplyr::filter(adj.P.Val < 0.05)

data.frame(
  Comparison = c("TaGL", "T1GH", "T2-4"),
  Significant_genes = c(
    nrow(sig_gene_TaGL),
    nrow(sig_gene_T1GH),
    nrow(sig_gene_T2_4)
  )
)

#Prepare the background gene universe for Enrichment Analysis
gene_universe <- gene_annotation_unique %>%
  dplyr::filter(
    !is.na(ENTREZID),
    ENTREZID != ""
  ) %>%
  dplyr::distinct(ENTREZID) %>%
  dplyr::pull(ENTREZID)

length(gene_universe)

#Prepare the significant gene lists 
genes_TaGL <- sig_gene_TaGL$ENTREZID
genes_T1GH <- sig_gene_T1GH$ENTREZID
genes_T2_4 <- sig_gene_T2_4$ENTREZID

genes_TaGL <- unique(na.omit(genes_TaGL))
genes_T1GH <- unique(na.omit(genes_T1GH))
genes_T2_4 <- unique(na.omit(genes_T2_4))

############ GO Biological Process enrichment analysis ##############
ego_TaGL <- enrichGO(
  gene = genes_TaGL,
  universe = gene_universe,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

ego_T1GH <- enrichGO(
  gene = genes_T1GH,
  universe = gene_universe,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

ego_T2_4 <- enrichGO(
  gene = genes_T2_4,
  universe = gene_universe,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

#Check the number of enriched biological processes
c(
  TaGL = nrow(as.data.frame(ego_TaGL)),
  T1GH = nrow(as.data.frame(ego_T1GH)),
  T2_4 = nrow(as.data.frame(ego_T2_4))
)

#Visualize the enriched GO Biological Processes for each comparison
dotplot(
  ego_TaGL,
  showCategory = 15,
  title = "GO Biological Processes - TaGL vs Control"
)

dotplot(
  ego_T1GH,
  showCategory = 15,
  title = "GO Biological Processes - T1GH vs Control"
)

dotplot(
  ego_T2_4,
  showCategory = 15,
  title = "GO Biological Processes - T2-4 vs Control"
)

#Save 
# GO BP: TaGL
png("GO_BP_TaGL.png", width = 3000, height = 2400, res = 300)
print(dotplot(ego_TaGL, showCategory = 15))
dev.off()


# GO BP: T1GH
png("GO_BP_T1GH.png", width = 3000, height = 2400, res = 300)
print(dotplot(ego_T1GH, showCategory = 15))
dev.off()


# GO BP: T2-4
png("GO_BP_T2_4.png", width = 3000, height = 2400, res = 300)
print(dotplot(ego_T2_4, showCategory = 15))
dev.off()

#Compare the enriched biological processes across the three comparisons
gene_list <- list(
  TaGL = genes_TaGL,
  T1GH = genes_T1GH,
  T2_4 = genes_T2_4
)

compare_go <- compareCluster(
  geneCluster = gene_list,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  universe = gene_universe
)

go_bp_comparison_plot <- dotplot(
  compare_go,
  showCategory = 15,
  title = "GO Biological Process Comparison",
  font.size = 8
)

#Save 
png(
  "GO_BP_Comparison.png",
  width = 3600,
  height = 2700,
  res = 300
)

print(go_bp_comparison_plot)

dev.off()

################ KEGG pathway enrichment analysis ###################
kegg_TaGL <- enrichKEGG(
  gene = genes_TaGL,
  universe = gene_universe,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

kegg_T1GH <- enrichKEGG(
  gene = genes_T1GH,
  universe = gene_universe,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

kegg_T2_4 <- enrichKEGG(
  gene = genes_T2_4,
  universe = gene_universe,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

c(
  TaGL = nrow(as.data.frame(kegg_TaGL)),
  T1GH = nrow(as.data.frame(kegg_T1GH)),
  T2_4 = nrow(as.data.frame(kegg_T2_4))
)

#Visualize 
dotplot(
  kegg_TaGL,
  showCategory = 15,
  title = "KEGG Pathways - TaGL vs Control",
  font.size = 9
)

dotplot(
  kegg_T1GH,
  showCategory = 15,
  title = "KEGG Pathways - T1GH vs Control",
  font.size = 9
)

dotplot(
  kegg_T2_4,
  showCategory = 15,
  title = "KEGG Pathways - T2-4 vs Control",
  font.size = 9
)

#Save
# KEGG: TaGL
png("KEGG_TaGL.png", width = 3000, height = 2400, res = 300)
print(dotplot(kegg_TaGL, showCategory = 15, font.size = 9))
dev.off()


# KEGG: T1GH
png("KEGG_T1GH.png", width = 3000, height = 2400, res = 300)
print(dotplot(kegg_T1GH, showCategory = 15, font.size = 9))
dev.off()


# KEGG: T2-4
png("KEGG_T2_4.png", width = 3000, height = 2400, res = 300)
print(dotplot(kegg_T2_4, showCategory = 15, font.size = 9))
dev.off()

#Compare KEGG pathways across the three conditions
compare_kegg <- compareCluster(
  geneCluster = gene_list,
  fun = "enrichKEGG",
  organism = "hsa",
  keyType = "ncbi-geneid",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  universe = gene_universe
)

#Visualize 
dotplot(
  compare_kegg,
  showCategory = 15,
  title = "KEGG Pathway Comparison",
  font.size = 9
)

#Save 
png(
  "KEGG_Pathway_Comparison.png",
  width = 3600,
  height = 2700,
  res = 300
)

print(
  dotplot(
    compare_kegg,
    showCategory = 15,
    title = "KEGG Pathway Comparison",
    font.size = 9
  )
)

dev.off()

############ Reactome pathway enrichment ##############
reactome_TaGL <- enrichPathway(
  gene = genes_TaGL,
  universe = gene_universe,
  organism = "human",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

reactome_T1GH <- enrichPathway(
  gene = genes_T1GH,
  universe = gene_universe,
  organism = "human",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

reactome_T2_4 <- enrichPathway(
  gene = genes_T2_4,
  universe = gene_universe,
  organism = "human",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

#Check the number of enriched pathways 
c(
  TaGL = nrow(as.data.frame(reactome_TaGL)),
  T1GH = nrow(as.data.frame(reactome_T1GH)),
  T2_4 = nrow(as.data.frame(reactome_T2_4))
)

#Visualize 
dotplot(
  reactome_TaGL,
  showCategory = 15,
  title = "Reactome Pathways - TaGL vs Control",
  font.size = 9
)

dotplot(
  reactome_T1GH,
  showCategory = 15,
  title = "Reactome Pathways - T1GH vs Control",
  font.size = 9
)

dotplot(
  reactome_T2_4,
  showCategory = 15,
  title = "Reactome Pathways - T2-4 vs Control",
  font.size = 9
)

#Save
# Reactome: TaGL
png("Reactome_TaGL.png", width = 3000, height = 2400, res = 300)
print(dotplot(reactome_TaGL, showCategory = 15, font.size = 9))
dev.off()


# Reactome: T1GH
png("Reactome_T1GH.png", width = 3000, height = 2400, res = 300)
print(dotplot(reactome_T1GH, showCategory = 15, font.size = 9))
dev.off()


# Reactome: T2-4
png("Reactome_T2_4.png", width = 3000, height = 2400, res = 300)
print(dotplot(reactome_T2_4, showCategory = 15, font.size = 9))
dev.off()

# Compare Reactome pathways across the three conditions
#Perform the combined analysis
compare_reactome <- compareCluster(
  geneCluster = gene_list,
  fun = "enrichPathway",
  organism = "human",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  universe = gene_universe
)

#Visualize the results
dotplot(
  compare_reactome,
  showCategory = 15,
  title = "Reactome Pathway Comparison",
  font.size = 7
)

#Save 
png(
  "Reactome_Pathway_Comparison.png",
  width = 3600,
  height = 2700,
  res = 300
)

print(
  dotplot(
    compare_reactome,
    showCategory = 15,
    title = "Reactome Pathway Comparison",
    font.size = 9
  )
)

dev.off()

############ Integrated pathway analysis ##############
#Export all enrichment results
extract_enrichment <- function(x, comparison, database) {
  
  df <- base::as.data.frame(x)
  
  df$Comparison <- comparison
  df$Database <- database
  
  return(df)
}

all_enrichment <- bind_rows(
  extract_enrichment(ego_TaGL, "TaGL", "GO-BP"),
  extract_enrichment(ego_T1GH, "T1GH", "GO-BP"),
  extract_enrichment(ego_T2_4, "T2-4", "GO-BP"),
  
  extract_enrichment(kegg_TaGL, "TaGL", "KEGG"),
  extract_enrichment(kegg_T1GH, "T1GH", "KEGG"),
  extract_enrichment(kegg_T2_4, "T2-4", "KEGG"),
  
  extract_enrichment(reactome_TaGL, "TaGL", "Reactome"),
  extract_enrichment(reactome_T1GH, "T1GH", "Reactome"),
  extract_enrichment(reactome_T2_4, "T2-4", "Reactome")
)

#Check the combined table
dim(all_enrichment)

table(all_enrichment$Database, all_enrichment$Comparison)

#Save 
write.csv(
  all_enrichment,
  "Integrated_Enrichment_Results.csv",
  row.names = FALSE
)

###################### Pathways shared across the three comparisons #######################
#Identify shared pathways
shared_pathways <- all_enrichment %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::group_by(Database, ID, Description) %>%
  dplyr::summarise(
    Number_of_comparisons = dplyr::n_distinct(Comparison),
    Comparisons = paste(
      sort(unique(Comparison)),
      collapse = ", "
    ),
    .groups = "drop"
  ) %>%
  dplyr::filter(Number_of_comparisons >= 2) %>%
  dplyr::arrange(Database, dplyr::desc(Number_of_comparisons))

#Inspect the results
head(shared_pathways, 20)

table(
  shared_pathways$Database,
  shared_pathways$Number_of_comparisons
)

#Separate upregulated and downregulated pathway enrichment
up_gene_list <- list(
  TaGL = sig_gene_TaGL %>%
    dplyr::filter(logFC > 0) %>%
    dplyr::pull(ENTREZID) %>%
    unique(),
  
  T1GH = sig_gene_T1GH %>%
    dplyr::filter(logFC > 0) %>%
    dplyr::pull(ENTREZID) %>%
    unique(),
  
  T2_4 = sig_gene_T2_4 %>%
    dplyr::filter(logFC > 0) %>%
    dplyr::pull(ENTREZID) %>%
    unique()
)

down_gene_list <- list(
  TaGL = sig_gene_TaGL %>%
    dplyr::filter(logFC < 0) %>%
    dplyr::pull(ENTREZID) %>%
    unique(),
  
  T1GH = sig_gene_T1GH %>%
    dplyr::filter(logFC < 0) %>%
    dplyr::pull(ENTREZID) %>%
    unique(),
  
  T2_4 = sig_gene_T2_4 %>%
    dplyr::filter(logFC < 0) %>%
    dplyr::pull(ENTREZID) %>%
    unique()
)

#Check the number of genes in each list
sapply(up_gene_list, length)
sapply(down_gene_list, length)


###################### Direction-specific GO enrichment #######################
# Upregulated genes
go_up <- clusterProfiler::compareCluster(
  geneCluster = up_gene_list,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  universe = gene_universe,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# Downregulated genes
go_down <- clusterProfiler::compareCluster(
  geneCluster = down_gene_list,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  universe = gene_universe,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

#Check the results 
nrow(as.data.frame(go_up))
nrow(as.data.frame(go_down))

#Compare upregulated and downregulated pathways
#Plot upregulated GO terms
up_go_plot <- enrichplot::dotplot(
  go_up,
  showCategory = 15,
  font.size = 10
) +
  ggplot2::ggtitle("GO Biological Process – Upregulated Genes") +
  ggplot2::theme_bw()

up_go_plot

#Plot downregulated GO terms
down_go_plot <- enrichplot::dotplot(
  go_down,
  showCategory = 15,
  font.size = 10
) +
  ggplot2::ggtitle("GO Biological Process – Downregulated Genes") +
  ggplot2::theme_bw()

down_go_plot

#Save
# GO upregulated comparison
png(
  "GO_Upregulated_Comparison.png",
  width = 3600,
  height = 2700,
  res = 300
)

print(up_go_plot)
dev.off()


# GO downregulated comparison
png(
  "GO_Downregulated_Comparison.png",
  width = 3600,
  height = 2700,
  res = 300
)

print(down_go_plot)
dev.off()

###################### Direction-specific KEGG enrichment #######################
#Upregulated genes
kegg_up <- clusterProfiler::compareCluster(
  geneCluster = up_gene_list,
  fun = "enrichKEGG",
  organism = "hsa",
  keyType = "ncbi-geneid",
  universe = gene_universe,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

#Downregulated genes
kegg_down <- clusterProfiler::compareCluster(
  geneCluster = down_gene_list,
  fun = "enrichKEGG",
  organism = "hsa",
  keyType = "ncbi-geneid",
  universe = gene_universe,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

#Check results 
nrow(as.data.frame(kegg_up))
nrow(as.data.frame(kegg_down))

#Visualize 
#Upregulated genes
kegg_up_plot <- enrichplot::dotplot(
  kegg_up,
  showCategory = 15,
  font.size = 10
) +
  ggplot2::ggtitle("KEGG Pathways – Upregulated Genes") +
  ggplot2::theme_bw()

kegg_up_plot

#Downregulated genes
kegg_down_plot <- enrichplot::dotplot(
  kegg_down,
  showCategory = 15,
  font.size = 10
) +
  ggplot2::ggtitle("KEGG Pathways – Downregulated Genes") +
  ggplot2::theme_bw()

kegg_down_plot

#Save
# KEGG upregulated comparison
png(
  "KEGG_Upregulated_Comparison.png",
  width = 3600,
  height = 2700,
  res = 300
)

print(kegg_up_plot)
dev.off()


# KEGG downregulated comparison
png(
  "KEGG_Downregulated_Comparison.png",
  width = 3600,
  height = 2700,
  res = 300
)

print(kegg_down_plot)
dev.off()

###################### Direction-specific Reactome enrichment #######################
#Upregulated genes
reactome_up <- clusterProfiler::compareCluster(
  geneCluster = up_gene_list,
  fun = "enrichPathway",
  organism = "human",
  universe = gene_universe,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

#Downregulated genes
reactome_down <- clusterProfiler::compareCluster(
  geneCluster = down_gene_list,
  fun = "enrichPathway",
  organism = "human",
  universe = gene_universe,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

#Check the results
nrow(as.data.frame(reactome_up))
nrow(as.data.frame(reactome_down))

#Visualize
#Upregulated genes
reactome_up_plot <- enrichplot::dotplot(
  reactome_up,
  showCategory = 15,
  font.size = 10
) +
  ggplot2::ggtitle("Reactome Pathways – Upregulated Genes") +
  ggplot2::theme_bw()

reactome_up_plot

#Downregulated genes
reactome_down_plot <- enrichplot::dotplot(
  reactome_down,
  showCategory = 15,
  font.size = 10
) +
  ggplot2::ggtitle("Reactome Pathways – Downregulated Genes") +
  ggplot2::theme_bw()

reactome_down_plot

#Save
# Reactome upregulated comparison
png(
  "Reactome_Upregulated_Comparison.png",
  width = 3600,
  height = 2700,
  res = 300
)

print(reactome_up_plot)
dev.off()


# Reactome downregulated comparison
png(
  "Reactome_Downregulated_Comparison.png",
  width = 3600,
  height = 2700,
  res = 300
)

print(reactome_down_plot)
dev.off()

###################### Prepare the heatmap expression matrix ########################
#Top 50 significant genes from each comparison, selected by adjusted p-value
# Select top 50 genes from each comparison
top_genes_TaGL <- sig_gene_TaGL %>%
  dplyr::arrange(adj.P.Val) %>%
  dplyr::slice_head(n = 50) %>%
  dplyr::pull(PROBEID)

top_genes_T1GH <- sig_gene_T1GH %>%
  dplyr::arrange(adj.P.Val) %>%
  dplyr::slice_head(n = 50) %>%
  dplyr::pull(PROBEID)

top_genes_T2_4 <- sig_gene_T2_4 %>%
  dplyr::arrange(adj.P.Val) %>%
  dplyr::slice_head(n = 50) %>%
  dplyr::pull(PROBEID)

#Combine them
top_genes_all <- unique(c(
  top_genes_TaGL,
  top_genes_T1GH,
  top_genes_T2_4
))

length(top_genes_all)

#Extract and standardize expression data
#Create the expression matrix
heatmap_matrix <- exp_norm[
  rownames(exp_norm) %in% top_genes_all,
  ,
  drop = FALSE
]

# Ensure samples are in metadata order
heatmap_matrix <- heatmap_matrix[
  ,
  pheno_data$Sample.Id,
  drop = FALSE
]

#Convert expression values to row-wise Z-scores
heatmap_matrix_scaled <- t(scale(t(heatmap_matrix)))

#Remove rows with undefined values
heatmap_matrix_scaled <- heatmap_matrix_scaled[
  complete.cases(heatmap_matrix_scaled),
  ,
  drop = FALSE
]

dim(heatmap_matrix_scaled)

#Prepare gene labels
heatmap_annotation <- gene_annotation_unique[
  match(rownames(heatmap_matrix_scaled),
        gene_annotation_unique$PROBEID),
]

gene_labels <- ifelse(
  is.na(heatmap_annotation$SYMBOL),
  rownames(heatmap_matrix_scaled),
  paste0(
    heatmap_annotation$SYMBOL,
    " (",
    rownames(heatmap_matrix_scaled),
    ")"
  )
)

rownames(heatmap_matrix_scaled) <- make.unique(gene_labels)

#Prepare sample annotations
sample_annotation <- data.frame(
  Sample.Type = pheno_data$Sample.Type
)

rownames(sample_annotation) <- pheno_data$Sample.Id

top_ha <- ComplexHeatmap::HeatmapAnnotation(
  df = sample_annotation
)

#Heatmap
heatmap_plot <- ComplexHeatmap::Heatmap(
  heatmap_matrix_scaled,
  name = "Z-score",
  top_annotation = top_ha,
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_title = "Top 93 Differentially Expressed Probe Sets",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  clustering_method_rows = "complete",
  clustering_method_columns = "complete"
)

heatmap_plot

# Save heatmap as high-resolution PNG (300 dpi)
png(
  "Top_93_DEG_Heatmap.png",
  width = 3000,
  height = 3600,
  res = 300
)

ComplexHeatmap::draw(heatmap_plot)

dev.off()

################## Identify common significant genes #####################
#Extract shared gene symbols
shared_genes <- Reduce(
  intersect,
  list(
    sig_gene_TaGL$SYMBOL,
    sig_gene_T1GH$SYMBOL,
    sig_gene_T2_4$SYMBOL
  )
)

shared_genes <- shared_genes[!is.na(shared_genes)]

length(shared_genes)

#Display the shared genes
head(shared_genes, 20)

#Create a table containing their statistics in all three comparisons
shared_gene_table <- dplyr::bind_rows(
  sig_gene_TaGL %>%
    dplyr::filter(SYMBOL %in% shared_genes) %>%
    dplyr::mutate(Comparison = "TaGL"),
  
  sig_gene_T1GH %>%
    dplyr::filter(SYMBOL %in% shared_genes) %>%
    dplyr::mutate(Comparison = "T1GH"),
  
  sig_gene_T2_4 %>%
    dplyr::filter(SYMBOL %in% shared_genes) %>%
    dplyr::mutate(Comparison = "T2-4")
)

#Inspect the table
dim(shared_gene_table)

head(shared_gene_table)

#Classify shared genes by regulation pattern
shared_gene_direction <- shared_gene_table %>%
  dplyr::select(SYMBOL, Comparison, logFC) %>%
  dplyr::mutate(
    Direction = dplyr::if_else(
      logFC > 0,
      "Upregulated",
      "Downregulated"
    )
  ) %>%
  dplyr::select(SYMBOL, Comparison, Direction) %>%
  tidyr::pivot_wider(
    names_from = Comparison,
    values_from = Direction
  )

shared_gene_direction <- shared_gene_direction %>%
  dplyr::mutate(
    Regulation_Pattern = dplyr::case_when(
      TaGL == "Upregulated" &
        T1GH == "Upregulated" &
        `T2-4` == "Upregulated" ~ "Consistently upregulated",
      
      TaGL == "Downregulated" &
        T1GH == "Downregulated" &
        `T2-4` == "Downregulated" ~ "Consistently downregulated",
      
      TRUE ~ "Oppositely regulated"
    )
  )

#Count each category
table(shared_gene_direction$Regulation_Pattern)

#Calculate percentages
prop.table(
  table(shared_gene_direction$Regulation_Pattern)
) * 100

#Identify the three oppositely regulated genes
opposite_genes <- shared_gene_direction %>%
  dplyr::filter(
    Regulation_Pattern == "Oppositely regulated"
  )

opposite_genes

#Identify their fold changes and adjusted p-values
shared_gene_table %>%
  dplyr::filter(SYMBOL %in% opposite_genes$SYMBOL) %>%
  dplyr::select(
    SYMBOL, Comparison, logFC, adj.P.Val
  ) %>%
  dplyr::arrange(SYMBOL, Comparison)

#Prioritize consistently regulated candidate genes
candidate_gene_summary <- shared_gene_table %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::summarise(
    Mean_logFC = mean(logFC),
    Mean_abs_logFC = mean(abs(logFC)),
    Max_adj_P = max(adj.P.Val),
    Mean_adj_P = mean(adj.P.Val),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    shared_gene_direction %>%
      dplyr::select(SYMBOL, Regulation_Pattern),
    by = "SYMBOL"
  ) %>%
  dplyr::arrange(Mean_adj_P)

#Separate the candidates
upregulated_candidates <- candidate_gene_summary %>%
  dplyr::filter(
    Regulation_Pattern == "Consistently upregulated"
  )

downregulated_candidates <- candidate_gene_summary %>%
  dplyr::filter(
    Regulation_Pattern == "Consistently downregulated"
  )

#Inspect the top 20 candidates in each group
head(upregulated_candidates, 20)

head(downregulated_candidates, 20)

##Select the strongest candidate genes
#Rank upregulated candidates
top_upregulated <- upregulated_candidates %>%
  dplyr::arrange(
    Max_adj_P,
    dplyr::desc(Mean_abs_logFC)
  ) %>%
  dplyr::slice_head(n = 20)

#Rank downregulated condidates
top_downregulated <- downregulated_candidates %>%
  dplyr::arrange(
    Max_adj_P,
    dplyr::desc(Mean_abs_logFC)
  ) %>%
  dplyr::slice_head(n = 20)

#Combine the results 
top_candidate_genes <- dplyr::bind_rows(
  top_upregulated,
  top_downregulated
)

top_candidate_genes

#Save the candidate table
write.csv(
  top_candidate_genes,
  "Top_40_Candidate_Genes.csv",
  row.names = FALSE
)

################ Gene Set Enrichment Analysis (GSEA) ##################
# Create a function to prepare ranked gene lists
prepare_gsea_ranks <- function(deg_table, annotation) {
  
  ranked_data <- deg_table %>%
    dplyr::left_join(
      annotation %>%
        dplyr::select(PROBEID, SYMBOL),
      by = "PROBEID"
    ) %>%
    dplyr::filter(
      !is.na(SYMBOL),
      SYMBOL != "",
      !is.na(t)
    ) %>%
    dplyr::arrange(dplyr::desc(abs(t))) %>%
    dplyr::distinct(SYMBOL, .keep_all = TRUE)
  
  gene_ranks <- ranked_data$t
  names(gene_ranks) <- ranked_data$SYMBOL
  
  sort(gene_ranks, decreasing = TRUE)
}

#Create the three ranked lists
ranks_TaGL <- prepare_gsea_ranks(
  deg_TaGL,
  gene_annotation_unique
)

ranks_T1GH <- prepare_gsea_ranks(
  deg_T1GH,
  gene_annotation_unique
)

ranks_T2_4 <- prepare_gsea_ranks(
  deg_T2_4,
  gene_annotation_unique
)

#Check the results 
length(ranks_TaGL)
length(ranks_T1GH)
length(ranks_T2_4)

head(ranks_TaGL)

#Retrieve Hallmark gene sets
hallmark <- msigdbr::msigdbr(
  species = "Homo sapiens",
  collection = "H"
)

head(hallmark)

#Prepare the Hallmark gene set database
hallmark_t2g <- hallmark %>%
  dplyr::select(gs_name, gene_symbol) %>%
  dplyr::distinct()

#Run GSEA for each comparison
gsea_TaGL <- clusterProfiler::GSEA(
  geneList = ranks_TaGL,
  TERM2GENE = hallmark_t2g,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = FALSE,
  seed = TRUE
)

gsea_T1GH <- clusterProfiler::GSEA(
  geneList = ranks_T1GH,
  TERM2GENE = hallmark_t2g,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = FALSE,
  seed = TRUE
)

gsea_T2_4 <- clusterProfiler::GSEA(
  geneList = ranks_T2_4,
  TERM2GENE = hallmark_t2g,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = FALSE,
  seed = TRUE
)

#Check the results
nrow(as.data.frame(gsea_TaGL))
nrow(as.data.frame(gsea_T1GH))
nrow(as.data.frame(gsea_T2_4))

#Compare GSEA results across conditions
gsea_TaGL_df <- as.data.frame(gsea_TaGL) %>%
  dplyr::mutate(Comparison = "TaGL")

gsea_T1GH_df <- as.data.frame(gsea_T1GH) %>%
  dplyr::mutate(Comparison = "T1GH")

gsea_T2_4_df <- as.data.frame(gsea_T2_4) %>%
  dplyr::mutate(Comparison = "T2-4")

all_gsea <- dplyr::bind_rows(
  gsea_TaGL_df,
  gsea_T1GH_df,
  gsea_T2_4_df
)

#Identify pathways shared across comparisons
shared_gsea_pathways <- all_gsea %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::group_by(ID, Description) %>%
  dplyr::summarise(
    Number_of_comparisons = dplyr::n_distinct(Comparison),
    Comparisons = paste(
      sort(unique(Comparison)),
      collapse = ", "
    ),
    .groups = "drop"
  ) %>%
  dplyr::filter(Number_of_comparisons >= 2) %>%
  dplyr::arrange(
    dplyr::desc(Number_of_comparisons)
  )

nrow(shared_gsea_pathways)

head(shared_gsea_pathways, 20)

#Classify shared pathways by enrichment direction
#Extract NES values for shared pathways
shared_gsea_details <- all_gsea %>%
  dplyr::filter(
    ID %in% shared_gsea_pathways$ID,
    p.adjust < 0.05
  ) %>%
  dplyr::select(
    ID, Description, Comparison, NES, p.adjust
  )

#Create a pathway comparison table
shared_gsea_wide <- shared_gsea_details %>%
  dplyr::select(ID, Description, Comparison, NES) %>%
  tidyr::pivot_wider(
    names_from = Comparison,
    values_from = NES
  )

#Classify pathway direction
shared_gsea_wide <- shared_gsea_wide %>%
  dplyr::mutate(
    Enrichment_Pattern = dplyr::case_when(
      TaGL > 0 & T1GH > 0 & `T2-4` > 0 ~
        "Consistently positively enriched",
      
      TaGL < 0 & T1GH < 0 & `T2-4` < 0 ~
        "Consistently negatively enriched",
      
      TRUE ~ "Oppositely enriched or missing"
    )
  )

#Count the patterns
table(shared_gsea_wide$Enrichment_Pattern)

#Inspect the 10 ambiguous pathways
shared_gsea_wide %>%
  dplyr::filter(
    Enrichment_Pattern == "Oppositely enriched or missing"
  ) %>%
  dplyr::select(
    Description, TaGL, T1GH, `T2-4`
  )

#GSEA visualization
#TaGL vs control
gsea_plot_TaGL <- enrichplot::dotplot(
  gsea_TaGL,
  showCategory = 20,
  font.size = 10
) +
  ggplot2::ggtitle("Hallmark GSEA – TaGL vs Control") +
  ggplot2::theme_bw()

gsea_plot_TaGL

#T1GH vs control
gsea_plot_T1GH <- enrichplot::dotplot(
  gsea_T1GH,
  showCategory = 20,
  font.size = 10
) +
  ggplot2::ggtitle("Hallmark GSEA – T1GH vs Control") +
  ggplot2::theme_bw()

gsea_plot_T1GH

#T2-4 vs control
gsea_plot_T2_4 <- enrichplot::dotplot(
  gsea_T2_4,
  showCategory = 20,
  font.size = 10
) +
  ggplot2::ggtitle("Hallmark GSEA – T2-4 vs Control") +
  ggplot2::theme_bw()

gsea_plot_T2_4

#Save
# GSEA: TaGL
png("GSEA_TaGL.png", width = 3300, height = 2700, res = 300)
print(gsea_plot_TaGL)
dev.off()


# GSEA: T1GH
png("GSEA_T1GH.png", width = 3300, height = 2700, res = 300)
print(gsea_plot_T1GH)
dev.off()


# GSEA: T2-4
png("GSEA_T2_4.png", width = 3300, height = 2700, res = 300)
print(gsea_plot_T2_4)
dev.off()

#Export GSEA results
write.csv(
  all_gsea,
  "All_Hallmark_GSEA_Results.csv",
  row.names = FALSE
)

write.csv(
  shared_gsea_wide,
  "Shared_Hallmark_GSEA_Pathways.csv",
  row.names = FALSE
)

#Retain all GSEA pathways
gsea_TaGL <- clusterProfiler::GSEA(
  geneList = ranks_TaGL,
  TERM2GENE = hallmark_t2g,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1,
  verbose = FALSE,
  seed = TRUE
)

gsea_T1GH <- clusterProfiler::GSEA(
  geneList = ranks_T1GH,
  TERM2GENE = hallmark_t2g,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1,
  verbose = FALSE,
  seed = TRUE
)

gsea_T2_4 <- clusterProfiler::GSEA(
  geneList = ranks_T2_4,
  TERM2GENE = hallmark_t2g,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1,
  verbose = FALSE,
  seed = TRUE
)

c(
  TaGL = nrow(as.data.frame(gsea_TaGL)),
  T1GH = nrow(as.data.frame(gsea_T1GH)),
  T2_4 = nrow(as.data.frame(gsea_T2_4))
)

all_gsea <- dplyr::bind_rows(
  as.data.frame(gsea_TaGL) %>% dplyr::mutate(Comparison = "TaGL"),
  as.data.frame(gsea_T1GH) %>% dplyr::mutate(Comparison = "T1GH"),
  as.data.frame(gsea_T2_4) %>% dplyr::mutate(Comparison = "T2-4")
)

#Identify Hallmark pathways that are statistically significant
shared_gsea_pathways <- all_gsea %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::group_by(ID, Description) %>%
  dplyr::summarise(
    Significant_Comparisons = dplyr::n_distinct(Comparison),
    .groups = "drop"
  ) %>%
  dplyr::filter(Significant_Comparisons >= 2) %>%
  dplyr::arrange(dplyr::desc(Significant_Comparisons))

nrow(shared_gsea_pathways)

table(shared_gsea_pathways$Significant_Comparisons)

#Compare the Normalized Enrichment Scores (NES) of the shared pathways
shared_gsea_details <- all_gsea %>%
  dplyr::filter(ID %in% shared_gsea_pathways$ID) %>%
  dplyr::select(ID, Description, Comparison, NES, p.adjust)

shared_gsea_wide <- shared_gsea_details %>%
  tidyr::pivot_wider(
    names_from = Comparison,
    values_from = c(NES, p.adjust)
  )

head(shared_gsea_wide)

#Classify each shared pathway according to its enrichment direction
shared_gsea_wide <- shared_gsea_wide %>%
  dplyr::mutate(
    Enrichment_Pattern = dplyr::case_when(
      NES_TaGL > 0 & NES_T1GH > 0 & `NES_T2-4` > 0 ~
        "Consistently positively enriched",
      
      NES_TaGL < 0 & NES_T1GH < 0 & `NES_T2-4` < 0 ~
        "Consistently negatively enriched",
      
      TRUE ~ "Oppositely enriched"
    )
  )

table(shared_gsea_wide$Enrichment_Pattern)

#Inspect the pathways in each enrichment category
shared_gsea_wide %>%
  dplyr::select(
    Description,
    NES_TaGL,
    NES_T1GH,
    `NES_T2-4`,
    Enrichment_Pattern
  ) %>%
  dplyr::arrange(Enrichment_Pattern, Description)

#Save the updated GSEA results and shared pathway comparison
write.csv(
  all_gsea,
  "All_Hallmark_GSEA_Results_Complete.csv",
  row.names = FALSE
)

write.csv(
  shared_gsea_pathways,
  "Shared_Hallmark_GSEA_Pathways.csv",
  row.names = FALSE
)

write.csv(
  shared_gsea_wide,
  "Shared_Hallmark_GSEA_Direction_Comparison.csv",
  row.names = FALSE
)

##################### Validate the differential-expression design matrix ######################
# Check sample group distribution
table(pheno_data$Sample.Type)

# Recreate the design matrix
design_check <- model.matrix(
  ~ 0 + Sample.Type,
  data = pheno_data
)

colnames(design_check)

# Inspect the design matrix
design_check

#Refit limma with verified contrasts
design <- design_check

colnames(design) <- c(
  "control",
  "TaGL",
  "T1GH",
  "T2_4"
)

fit_check <- limma::lmFit(exp_norm, design)

contrasts_check <- limma::makeContrasts(
  TaGL_vs_control = TaGL - control,
  T1GH_vs_control = T1GH - control,
  T2_4_vs_control = T2_4 - control,
  levels = design
)

contrasts_check

#Verify the contrast coefficients and refit the empirical Bayes model
fit2_check <- limma::contrasts.fit(
  fit_check,
  contrasts_check
)

fit2_check <- limma::eBayes(fit2_check)

# Check the fitted model
dim(fit2_check$coefficients)

colnames(fit2_check$coefficients)

summary(fit2_check)

#Compare the original and refitted differential-expression results
# Extract updated DEG results
deg_TaGL_check <- limma::topTable(
  fit2_check,
  coef = "TaGL_vs_control",
  number = Inf,
  adjust.method = "BH"
)

deg_T1GH_check <- limma::topTable(
  fit2_check,
  coef = "T1GH_vs_control",
  number = Inf,
  adjust.method = "BH"
)

deg_T2_4_check <- limma::topTable(
  fit2_check,
  coef = "T2_4_vs_control",
  number = Inf,
  adjust.method = "BH"
)

# Compare dimensions with original DEG tables
rbind(
  Original = c(
    TaGL = nrow(deg_TaGL),
    T1GH = nrow(deg_T1GH),
    T2_4 = nrow(deg_T2_4)
  ),
  Refit = c(
    TaGL = nrow(deg_TaGL_check),
    T1GH = nrow(deg_T1GH_check),
    T2_4 = nrow(deg_T2_4_check)
  )
)

#Compare logFC and adjusted p-values
# Compare original and refitted results
comparison_check <- data.frame(
  Comparison = c("TaGL", "T1GH", "T2-4"),
  
  Max_Abs_Difference_logFC = c(
    max(abs(deg_TaGL$logFC - deg_TaGL_check$logFC), na.rm = TRUE),
    max(abs(deg_T1GH$logFC - deg_T1GH_check$logFC), na.rm = TRUE),
    max(abs(deg_T2_4$logFC - deg_T2_4_check$logFC), na.rm = TRUE)
  ),
  
  Max_Abs_Difference_adjP = c(
    max(abs(deg_TaGL$adj.P.Val - deg_TaGL_check$adj.P.Val), na.rm = TRUE),
    max(abs(deg_T1GH$adj.P.Val - deg_T1GH_check$adj.P.Val), na.rm = TRUE),
    max(abs(deg_T2_4$adj.P.Val - deg_T2_4_check$adj.P.Val), na.rm = TRUE)
  )
)

comparison_check

#Validate sample matching
# Inspect sample IDs
colnames(exp_norm)

pheno_data$Sample.Id

# Check whether sample IDs match in the same order
identical(
  colnames(exp_norm),
  pheno_data$Sample.Id
)

#Validate probe-to-gene annotation
# Annotation coverage
annotation_summary <- data.frame(
  Total_Probes = nrow(exp_norm),
  Annotated_Probes = sum(!is.na(gene_annotation_unique$SYMBOL)),
  Unique_Genes = dplyr::n_distinct(gene_annotation_unique$SYMBOL),
  Unique_Entrez_IDs = dplyr::n_distinct(
    gene_annotation_unique$ENTREZID,
    na.rm = TRUE
  )
)

annotation_summary

# Number of probes per gene
probe_per_gene <- gene_annotation_unique %>%
  dplyr::filter(!is.na(SYMBOL)) %>%
  dplyr::count(SYMBOL, name = "Number_of_Probes")

summary(probe_per_gene$Number_of_Probes)

#Check annotation ambiguity
probe_mapping_check <- gene_annotation %>%
  dplyr::filter(!is.na(SYMBOL)) %>%
  dplyr::group_by(PROBEID) %>%
  dplyr::summarise(
    Number_of_Gene_Symbols = dplyr::n_distinct(SYMBOL),
    Gene_Symbols = paste(unique(SYMBOL), collapse = "; "),
    .groups = "drop"
  )

# Summary of probe mapping ambiguity
table(probe_mapping_check$Number_of_Gene_Symbols)

# Inspect probes mapping to multiple gene symbols
ambiguous_probes <- probe_mapping_check %>%
  dplyr::filter(Number_of_Gene_Symbols > 1) %>%
  dplyr::arrange(dplyr::desc(Number_of_Gene_Symbols))

head(ambiguous_probes, 10)

#Quantify the ambiguity
ambiguity_summary <- data.frame(
  Total_Annotated_Probes = nrow(probe_mapping_check),
  Unambiguous_Probes = sum(
    probe_mapping_check$Number_of_Gene_Symbols == 1
  ),
  Ambiguous_Probes = sum(
    probe_mapping_check$Number_of_Gene_Symbols > 1
  ),
  Ambiguous_Percentage = mean(
    probe_mapping_check$Number_of_Gene_Symbols > 1
  ) * 100
)

ambiguity_summary

head(ambiguous_probes, 10)

#Create the clean annotation
unambiguous_annotation <- probe_mapping_check %>%
  dplyr::filter(Number_of_Gene_Symbols == 1) %>%
  dplyr::select(PROBEID)

gene_annotation_clean <- gene_annotation_unique %>%
  dplyr::semi_join(
    unambiguous_annotation,
    by = "PROBEID"
  )

c(
  Annotated_Probes = nrow(gene_annotation_clean),
  Unique_Genes = dplyr::n_distinct(gene_annotation_clean$SYMBOL)
)

#Rebuild gene-level DEG tables
create_gene_table_clean <- function(deg_table, annotation) {
  
  deg_table %>%
    dplyr::inner_join(annotation, by = "PROBEID") %>%
    dplyr::arrange(adj.P.Val) %>%
    dplyr::distinct(SYMBOL, .keep_all = TRUE)
}

gene_TaGL_clean <- create_gene_table_clean(
  deg_TaGL,
  gene_annotation_clean
)

gene_T1GH_clean <- create_gene_table_clean(
  deg_T1GH,
  gene_annotation_clean
)

gene_T2_4_clean <- create_gene_table_clean(
  deg_T2_4,
  gene_annotation_clean
)

# Check dimensions
rbind(
  TaGL = c(
    Original = nrow(gene_TaGL),
    Clean = nrow(gene_TaGL_clean)
  ),
  T1GH = c(
    Original = nrow(gene_T1GH),
    Clean = nrow(gene_T1GH_clean)
  ),
  T2_4 = c(
    Original = nrow(gene_T2_4),
    Clean = nrow(gene_T2_4_clean)
  )
)

#Compare significant gene counts
sig_gene_TaGL_clean <- gene_TaGL_clean %>%
  dplyr::filter(adj.P.Val < 0.05)

sig_gene_T1GH_clean <- gene_T1GH_clean %>%
  dplyr::filter(adj.P.Val < 0.05)

sig_gene_T2_4_clean <- gene_T2_4_clean %>%
  dplyr::filter(adj.P.Val < 0.05)

significance_comparison <- data.frame(
  Comparison = c("TaGL", "T1GH", "T2-4"),
  
  Original_Significant = c(
    nrow(sig_gene_TaGL),
    nrow(sig_gene_T1GH),
    nrow(sig_gene_T2_4)
  ),
  
  Clean_Significant = c(
    nrow(sig_gene_TaGL_clean),
    nrow(sig_gene_T1GH_clean),
    nrow(sig_gene_T2_4_clean)
  )
) %>%
  dplyr::mutate(
    Difference = Clean_Significant - Original_Significant,
    Percentage_Retained = round(
      Clean_Significant / Original_Significant * 100, 2
    )
  )

significance_comparison

#Compare DEG identity
compare_gene_overlap <- function(original, clean, comparison) {
  
  original_genes <- original %>%
    dplyr::filter(adj.P.Val < 0.05) %>%
    dplyr::pull(SYMBOL)
  
  clean_genes <- clean %>%
    dplyr::filter(adj.P.Val < 0.05) %>%
    dplyr::pull(SYMBOL)
  
  common <- intersect(original_genes, clean_genes)
  
  data.frame(
    Comparison = comparison,
    Original_DEGs = length(original_genes),
    Clean_DEGs = length(clean_genes),
    Shared_DEGs = length(common),
    Jaccard_Index = round(
      length(common) / length(union(original_genes, clean_genes)),
      4
    ),
    Retention_Percentage = round(
      length(common) / length(original_genes) * 100, 2
    )
  )
}

gene_overlap_sensitivity <- dplyr::bind_rows(
  compare_gene_overlap(gene_TaGL, gene_TaGL_clean, "TaGL"),
  compare_gene_overlap(gene_T1GH, gene_T1GH_clean, "T1GH"),
  compare_gene_overlap(gene_T2_4, gene_T2_4_clean, "T2-4")
)

gene_overlap_sensitivity

#Compare fold changes
compare_logFC <- function(original, clean, comparison) {
  
  original_sig <- original %>%
    dplyr::filter(adj.P.Val < 0.05) %>%
    dplyr::select(SYMBOL, logFC) %>%
    dplyr::rename(logFC_original = logFC)
  
  clean_sig <- clean %>%
    dplyr::filter(adj.P.Val < 0.05) %>%
    dplyr::select(SYMBOL, logFC) %>%
    dplyr::rename(logFC_clean = logFC)
  
  merged <- dplyr::inner_join(
    original_sig,
    clean_sig,
    by = "SYMBOL"
  )
  
  data.frame(
    Comparison = comparison,
    Pearson_Correlation = cor(
      merged$logFC_original,
      merged$logFC_clean,
      method = "pearson"
    ),
    Spearman_Correlation = cor(
      merged$logFC_original,
      merged$logFC_clean,
      method = "spearman"
    ),
    Max_Abs_Difference = max(
      abs(merged$logFC_original - merged$logFC_clean)
    )
  )
}

logFC_sensitivity <- dplyr::bind_rows(
  compare_logFC(gene_TaGL, gene_TaGL_clean, "TaGL"),
  compare_logFC(gene_T1GH, gene_T1GH_clean, "T1GH"),
  compare_logFC(gene_T2_4, gene_T2_4_clean, "T2-4")
)

logFC_sensitivity

#Investigate the largest fold-change differences
find_largest_logFC_differences <- function(original, clean, comparison) {
  
  original_sig <- original %>%
    dplyr::filter(adj.P.Val < 0.05) %>%
    dplyr::select(SYMBOL, logFC) %>%
    dplyr::rename(logFC_original = logFC)
  
  clean_sig <- clean %>%
    dplyr::filter(adj.P.Val < 0.05) %>%
    dplyr::select(SYMBOL, logFC) %>%
    dplyr::rename(logFC_clean = logFC)
  
  dplyr::inner_join(original_sig, clean_sig, by = "SYMBOL") %>%
    dplyr::mutate(
      Difference = logFC_clean - logFC_original,
      Absolute_Difference = abs(Difference),
      Comparison = comparison
    ) %>%
    dplyr::arrange(dplyr::desc(Absolute_Difference)) %>%
    dplyr::slice_head(n = 10)
}

largest_logFC_differences <- dplyr::bind_rows(
  find_largest_logFC_differences(
    gene_TaGL, gene_TaGL_clean, "TaGL"
  ),
  find_largest_logFC_differences(
    gene_T1GH, gene_T1GH_clean, "T1GH"
  ),
  find_largest_logFC_differences(
    gene_T2_4, gene_T2_4_clean, "T2-4"
  )
)

largest_logFC_differences

#Investigate RHBDL2 probe mapping
# RHBDL2 probes in original annotation
gene_annotation_unique %>%
  dplyr::filter(SYMBOL == "RHBDL2")

# RHBDL2 probes in clean annotation
gene_annotation_clean %>%
  dplyr::filter(SYMBOL == "RHBDL2")

deg_TaGL %>%
  dplyr::filter(PROBEID %in% gene_annotation_clean$PROBEID[
    gene_annotation_clean$SYMBOL == "RHBDL2"
  ]) %>%
  dplyr::select(PROBEID, logFC, P.Value, adj.P.Val)

#Inspect all RHBDL2 probes
gene_annotation %>%
  dplyr::filter(PROBEID %in% c(
    "1554895_a_at",
    "1552502_s_at",
    "1554897_s_at"
  )) %>%
  dplyr::select(PROBEID, SYMBOL, GENENAME, ENTREZID)

gene_TaGL %>%
  dplyr::filter(SYMBOL == "RHBDL2") %>%
  dplyr::select(PROBEID, SYMBOL, logFC, P.Value, adj.P.Val)

#Visualize RHBDL2 expression across the four groups
rhbdl2_probes <- c(
  "219489_s_at",
  "1552502_s_at",
  "1554895_a_at",
  "1554897_s_at"
)

rhbdl2_expression <- as.data.frame(
  t(exp_norm[rhbdl2_probes, , drop = FALSE])
)

rhbdl2_expression$Sample.Id <- rownames(rhbdl2_expression)

rhbdl2_expression <- rhbdl2_expression %>%
  dplyr::left_join(
    pheno_data %>%
      dplyr::select(Sample.Id, Sample.Type),
    by = "Sample.Id"
  )

rhbdl2_long <- rhbdl2_expression %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(rhbdl2_probes),
    names_to = "Probe",
    values_to = "Expression"
  )

ggplot2::ggplot(
  rhbdl2_long,
  ggplot2::aes(x = Sample.Type, y = Expression, fill = Sample.Type)
) +
  ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  ggplot2::geom_jitter(width = 0.12, size = 2) +
  ggplot2::facet_wrap(~ Probe, scales = "free_y") +
  ggplot2::theme_bw() +
  ggplot2::labs(
    title = "RHBDL2 Probe-Level Expression",
    x = "Experimental Group",
    y = "Normalized Expression"
  )

#Quantify the RHBDL2 expression patterns
rhbdl2_group_means <- rhbdl2_long %>%
  dplyr::group_by(Probe, Sample.Type) %>%
  dplyr::summarise(
    Mean_Expression = mean(Expression),
    SD_Expression = sd(Expression),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = Sample.Type,
    values_from = c(Mean_Expression, SD_Expression)
  )

rhbdl2_group_means

#Verify the statistical significance of all four probes
rhbdl2_limma <- limma::topTable(
  fit2_check,
  coef = "TaGL_vs_control",
  number = Inf,
  adjust.method = "BH"
) %>%
  tibble::rownames_to_column("PROBEID") %>%
  dplyr::filter(PROBEID %in% rhbdl2_probes) %>%
  dplyr::select(PROBEID, logFC, AveExpr, P.Value, adj.P.Val)

rhbdl2_limma

#Check the probe annotation and transcript coverage
rhbdl2_annotation <- gene_annotation %>%
  dplyr::filter(PROBEID %in% rhbdl2_probes) %>%
  dplyr::distinct(PROBEID, SYMBOL, GENENAME, ENTREZID)

rhbdl2_annotation

#Confirm the ambiguous probe mapping
rhbdl2_mapping_check <- rhbdl2_annotation %>%
  dplyr::group_by(PROBEID) %>%
  dplyr::summarise(
    Number_of_Genes = dplyr::n_distinct(SYMBOL),
    Mapped_Genes = paste(unique(SYMBOL), collapse = "; "),
    .groups = "drop"
  )

rhbdl2_mapping_check

#Inspect the current ranking function
prepare_gsea_ranks

#Create an alternative GSEA ranking function that averages the moderated t-statistics across probes mapping to the same gene
prepare_gsea_ranks_mean <- function(deg_table, annotation) {
  
  gene_ranks <- deg_table %>%
    dplyr::inner_join(
      annotation %>% dplyr::select(PROBEID, SYMBOL),
      by = "PROBEID"
    ) %>%
    dplyr::filter(
      !is.na(SYMBOL),
      SYMBOL != "",
      !is.na(t)
    ) %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(
      Mean_t = mean(t),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(Mean_t))
  
  ranks <- gene_ranks$Mean_t
  names(ranks) <- gene_ranks$SYMBOL
  
  sort(ranks, decreasing = TRUE)
}

ranks_TaGL_mean <- prepare_gsea_ranks_mean(
  deg_TaGL, gene_annotation_clean
)

ranks_T1GH_mean <- prepare_gsea_ranks_mean(
  deg_T1GH, gene_annotation_clean
)

ranks_T2_4_mean <- prepare_gsea_ranks_mean(
  deg_T2_4, gene_annotation_clean
)

length(ranks_TaGL_mean)
length(ranks_T1GH_mean)
length(ranks_T2_4_mean)

c(
  TaGL = anyDuplicated(names(ranks_TaGL_mean)),
  T1GH = anyDuplicated(names(ranks_T1GH_mean)),
  T2_4 = anyDuplicated(names(ranks_T2_4_mean))
)

#Run sensitivity GSEA
gsea_TaGL_mean <- clusterProfiler::GSEA(
  geneList = ranks_TaGL_mean,
  TERM2GENE = hallmark_t2g,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  seed = TRUE,
  verbose = FALSE
)

gsea_T1GH_mean <- clusterProfiler::GSEA(
  geneList = ranks_T1GH_mean,
  TERM2GENE = hallmark_t2g,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  seed = TRUE,
  verbose = FALSE
)

gsea_T2_4_mean <- clusterProfiler::GSEA(
  geneList = ranks_T2_4_mean,
  TERM2GENE = hallmark_t2g,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  seed = TRUE,
  verbose = FALSE
)

nrow(as.data.frame(gsea_TaGL_mean))
nrow(as.data.frame(gsea_T1GH_mean))
nrow(as.data.frame(gsea_T2_4_mean))

c(
  TaGL = nrow(as.data.frame(gsea_TaGL_mean)),
  T1GH = nrow(as.data.frame(gsea_T1GH_mean)),
  T2_4 = nrow(as.data.frame(gsea_T2_4_mean))
)

#Combine the results into one table
gsea_mean_all <- dplyr::bind_rows(
  as.data.frame(gsea_TaGL_mean) %>%
    dplyr::mutate(Comparison = "TaGL"),
  as.data.frame(gsea_T1GH_mean) %>%
    dplyr::mutate(Comparison = "T1GH"),
  as.data.frame(gsea_T2_4_mean) %>%
    dplyr::mutate(Comparison = "T2-4")
)

gsea_mean_wide <- gsea_mean_all %>%
  dplyr::select(ID, Description, Comparison, NES, p.adjust) %>%
  tidyr::pivot_wider(
    names_from = Comparison,
    values_from = c(NES, p.adjust)
  )

dim(gsea_mean_wide)

#Compare original and mean-based NES values
gsea_nes_comparison <- data.frame(
  Comparison = c("TaGL", "T1GH", "T2-4"),
  Pearson_NES = c(
    cor(
      gsea_TaGL@result$NES,
      gsea_TaGL_mean@result$NES,
      use = "complete.obs"
    ),
    cor(
      gsea_T1GH@result$NES,
      gsea_T1GH_mean@result$NES,
      use = "complete.obs"
    ),
    cor(
      gsea_T2_4@result$NES,
      gsea_T2_4_mean@result$NES,
      use = "complete.obs"
    )
  )
)

gsea_nes_comparison

#Compare significant pathways
gsea_significance_comparison <- dplyr::bind_rows(
  data.frame(
    Comparison = "TaGL",
    Original_Significant = sum(gsea_TaGL@result$p.adjust < 0.05, na.rm = TRUE),
    Mean_Significant = sum(gsea_TaGL_mean@result$p.adjust < 0.05, na.rm = TRUE)
  ),
  data.frame(
    Comparison = "T1GH",
    Original_Significant = sum(gsea_T1GH@result$p.adjust < 0.05, na.rm = TRUE),
    Mean_Significant = sum(gsea_T1GH_mean@result$p.adjust < 0.05, na.rm = TRUE)
  ),
  data.frame(
    Comparison = "T2-4",
    Original_Significant = sum(gsea_T2_4@result$p.adjust < 0.05, na.rm = TRUE),
    Mean_Significant = sum(gsea_T2_4_mean@result$p.adjust < 0.05, na.rm = TRUE)
  )
)

gsea_significance_comparison

#Identify pathways significant in both methods
compare_significant_gsea <- function(original, mean_based, comparison) {
  
  original_sig <- original@result %>%
    dplyr::filter(p.adjust < 0.05) %>%
    dplyr::select(ID, Description, NES_original = NES)
  
  mean_sig <- mean_based@result %>%
    dplyr::filter(p.adjust < 0.05) %>%
    dplyr::select(ID, NES_mean = NES)
  
  comparison_table <- dplyr::inner_join(
    original_sig,
    mean_sig,
    by = "ID"
  ) %>%
    dplyr::mutate(
      Comparison = comparison,
      Direction_Agreement = sign(NES_original) == sign(NES_mean)
    )
  
  comparison_table
}

shared_sig_TaGL <- compare_significant_gsea(
  gsea_TaGL, gsea_TaGL_mean, "TaGL"
)

shared_sig_T1GH <- compare_significant_gsea(
  gsea_T1GH, gsea_T1GH_mean, "T1GH"
)

shared_sig_T2_4 <- compare_significant_gsea(
  gsea_T2_4, gsea_T2_4_mean, "T2-4"
)

gsea_overlap_summary <- dplyr::bind_rows(
  shared_sig_TaGL,
  shared_sig_T1GH,
  shared_sig_T2_4
) %>%
  dplyr::group_by(Comparison) %>%
  dplyr::summarise(
    Shared_Significant_Pathways = dplyr::n(),
    Same_Direction = sum(Direction_Agreement),
    Opposite_Direction = sum(!Direction_Agreement),
    .groups = "drop"
  )

gsea_overlap_summary

#Quantify pathway overlap
gsea_overlap_metrics <- function(original, mean_based, comparison) {
  
  original_sig <- original@result %>%
    dplyr::filter(p.adjust < 0.05) %>%
    dplyr::pull(ID)
  
  mean_sig <- mean_based@result %>%
    dplyr::filter(p.adjust < 0.05) %>%
    dplyr::pull(ID)
  
  shared <- intersect(original_sig, mean_sig)
  
  data.frame(
    Comparison = comparison,
    Original_Only = length(setdiff(original_sig, mean_sig)),
    Mean_Only = length(setdiff(mean_sig, original_sig)),
    Shared = length(shared),
    Jaccard_Index = length(shared) /
      length(union(original_sig, mean_sig))
  )
}

gsea_overlap_metrics_all <- dplyr::bind_rows(
  gsea_overlap_metrics(gsea_TaGL, gsea_TaGL_mean, "TaGL"),
  gsea_overlap_metrics(gsea_T1GH, gsea_T1GH_mean, "T1GH"),
  gsea_overlap_metrics(gsea_T2_4, gsea_T2_4_mean, "T2-4")
)

gsea_overlap_metrics_all

#Examine pathways unique to each ranking method
gsea_unique_pathways <- function(original, mean_based, comparison) {
  
  original_sig <- original@result %>%
    dplyr::filter(p.adjust < 0.05) %>%
    dplyr::select(ID, Description, NES, p.adjust)
  
  mean_sig <- mean_based@result %>%
    dplyr::filter(p.adjust < 0.05) %>%
    dplyr::select(ID, Description, NES, p.adjust)
  
  original_only <- original_sig %>%
    dplyr::anti_join(mean_sig, by = "ID") %>%
    dplyr::mutate(Method = "Original only")
  
  mean_only <- mean_sig %>%
    dplyr::anti_join(original_sig, by = "ID") %>%
    dplyr::mutate(Method = "Mean-based only")
  
  dplyr::bind_rows(original_only, mean_only) %>%
    dplyr::mutate(Comparison = comparison)
}

gsea_unique_all <- dplyr::bind_rows(
  gsea_unique_pathways(gsea_TaGL, gsea_TaGL_mean, "TaGL"),
  gsea_unique_pathways(gsea_T1GH, gsea_T1GH_mean, "T1GH"),
  gsea_unique_pathways(gsea_T2_4, gsea_T2_4_mean, "T2-4")
)

gsea_unique_all %>%
  dplyr::select(
    Comparison, Description, Method, NES, p.adjust
  )

