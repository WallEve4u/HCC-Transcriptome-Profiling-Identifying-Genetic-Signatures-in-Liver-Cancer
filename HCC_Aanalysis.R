################################################################################
# Project: Differential Gene Expression Analysis (HCC vs Normal)
# Author: Waleed Tariq
# Dataset: GSE112790
################################################################################

#--- 1. Loading Essential Libraries ---
library(GEOquery)
library(limma)
library(hgu133plus2.db)
library(ggplot2)
library(ggrepel)
library(pheatmap)

#--- 2. Data Acquisition & Group Classification ---
# Using locally downloaded file — most reliable method
gse <- getGEO(filename = "GSE112790_series_matrix.txt.gz", getGPL = FALSE)
# Note: filename= returns ExpressionSet directly, NOT a list — no [[1]] needed

metadata <- pData(gse)
expr     <- exprs(gse)

# Defining experimental groups based on sample source description
group <- ifelse(
  grepl("normal", metadata$source_name_ch1, ignore.case = TRUE),
  "Normal", "Tumor"
)
group <- factor(group)

# Sanity check — run this to confirm groups look correct before proceeding
print(table(group))

#--- 3. Statistical Analysis (Differential Expression Analysis) ---
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

contrast <- makeContrasts(Tumor_vs_Normal = Tumor - Normal, levels = design)

fit  <- lmFit(expr, design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

sigGenes <- topTable(fit2, number = Inf, p.value = 0.05, lfc = 1)
cat("Number of significant DEGs found:", nrow(sigGenes), "\n")

#--- 4. Biological Annotation (Probe ID → Gene Symbol) ---
my_probes <- rownames(sigGenes)
geneMap   <- select(hgu133plus2.db,
                    keys    = my_probes,
                    columns = "SYMBOL",
                    keytype = "PROBEID")

final_results <- merge(sigGenes, geneMap, by.x = "row.names", by.y = "PROBEID")
colnames(final_results)[1] <- "ProbeID"
final_results <- final_results[order(final_results$adj.P.Val), ]

#--- 5. Volcano Plot ---
final_results$diffexpressed <- "NO"
final_results$diffexpressed[final_results$logFC >  1 & final_results$adj.P.Val < 0.05] <- "UP"
final_results$diffexpressed[final_results$logFC < -1 & final_results$adj.P.Val < 0.05] <- "DOWN"

ggplot(data = final_results,
       aes(x = logFC, y = -log10(adj.P.Val), col = diffexpressed)) +
  geom_point(alpha = 0.4, size = 1.5) +
  theme_minimal() +
  scale_color_manual(values = c("DOWN" = "blue", "UP" = "red", "NO" = "grey")) +
  geom_vline(xintercept = c(-1, 1),    color = "black", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), color = "black", linetype = "dashed") +
  geom_text_repel(data = head(final_results, 10), aes(label = SYMBOL), size = 4) +
  labs(title = "Volcano Plot: Liver Cancer Analysis",
       x = "log2 Fold Change", y = "-log10 Adj P-value")

#--- 6. Clustered Heatmap (Top 50 DEGs) ---
final_filtred <- final_results[!is.na(final_results$SYMBOL), ]
final_filtred <- final_filtred[!duplicated(final_filtred$SYMBOL), ]

top_50_ids        <- final_filtred$ProbeID[1:50]
plot_matrix       <- expr[top_50_ids, ]
rownames(plot_matrix) <- final_filtred$SYMBOL[1:50]
plot_matrix       <- t(scale(t(plot_matrix)))  # Z-score normalisation

sample_info <- data.frame(Group = group)
rownames(sample_info) <- colnames(expr)
ann_color   <- list(Group = c(Normal = "blue", Tumor = "red"))

pheatmap(plot_matrix,
         annotation_col    = sample_info,
         annotation_colors = ann_color,
         main              = "Heatmap: Top 50 Differentially Expressed Genes",
         color             = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         show_colnames     = FALSE,
         fontsize_row      = 8)

#--- 7. Data Export ---
rownames(final_results) <- NULL
write.csv(final_results, "HCC_Significant_Genes.csv", row.names = FALSE)
cat("Results saved to HCC_Significant_Genes.csv\n")