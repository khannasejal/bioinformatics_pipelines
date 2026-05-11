rm(list = ls()) # Clear all variables from the environment
gc()

library(Seurat)
library(Signac)
library(ggplot2)
library(stringr)
library(dplyr)
library(tidyr)

#seurat_merged_WT = readRDS("seurat_merged_WT.RDS")
split_plot <- DimPlot(seurat_merged_WT, 
                      reduction = "umap_wnn", 
                      split.by = "orig.ident") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
seurat_merged_WT = readRDS("seurat_merged_MUT.RDS")
#sample_order <- c("E","H1_0h","H1_24h","H1_48h","H1_72h","H2","M")
#sample_order <- c("WT_0h","WT_24h","WT_48h","WT_72h")
sample_order <- c("MUT_0h","MUT_24h","MUT_48h","MUT_72h")   #while running for wild type samples, change the sample order accordingly
#sample_order <- c("E","H1","H2","M")

DefaultAssay(seurat_merged_WT) <- "RNA"
Idents(seurat_merged_WT) <- "orig.ident" #don't run this for the cluster comparisons
DimPlot(seurat_merged_WT, reduction = "umap_wnn", label = TRUE)

seurat_ordered <- seurat_merged_WT

Idents(seurat_ordered) <- factor(
  Idents(seurat_ordered),
  levels = sample_order
)

final_table <- FetchData(
  object = seurat_ordered, 
  vars = c("orig.ident", "Vim", "Cdh1"), 
  layer = "data"
)

final_table$Vim_Cdh1_diff <- final_table$Vim - final_table$Cdh1
colnames(final_table)[1] <- "Sample"
head(final_table)

write.table(
  x = final_table, 
  file = "wt_vim_cdh1_exp_time_points.txt", 
  row.names = TRUE,
  quote = FALSE,# Keeps the unique Cell Barcodes
  sep = "\t"
)


gene_names <- rownames(seurat_merged_WT[["RNA"]])

write.table(
  x = gene_names, 
  file = "mouse_gene_names.txt", 
  sep = "\t", 
  quote = FALSE, 
  row.names = FALSE, 
  col.names = FALSE
)

FeaturePlot(
  object = seurat_ordered, 
  features = "Vim", 
  reduction = "umap_wnn", # Customizing colors: grey (low) to blue (high)
  label = TRUE                   # Optional: keep cluster labels on the plot
)


mes_genes <- str_to_title(read.table("./genesets/mes_cell_line_signature.txt", stringsAsFactors = FALSE)$V1)
#mes_genes <- str_to_title(read.table("./genesets/chiara_mes_listC.txt", stringsAsFactors = FALSE)$V1)
#mes_genes <- str_to_title(read.table("./genesets/chiara_mes_listB.txt", stringsAsFactors = FALSE)$V1)
#mes_genes <- str_to_title(read.table("./genesets/chiara_mes_listA.txt", stringsAsFactors = FALSE)$V1)


epi_genes <- str_to_title(read.table("./genesets/epi_cell_line_signature_w_cdh1.txt", stringsAsFactors = FALSE)$V1)
#epi_genes <- str_to_title(read.table("./genesets/chiara_epi_listC.txt", stringsAsFactors = FALSE)$V1)
#epi_genes <- str_to_title(read.table("./genesets/chiara_epi_listB.txt", stringsAsFactors = FALSE)$V1)
#epi_genes <- str_to_title(read.table("./genesets/chiara_epi_listA.txt", stringsAsFactors = FALSE)$V1)


mes_valid <- mes_genes[mes_genes %in% rownames(seurat_ordered)]
epi_valid <- epi_genes[epi_genes %in% rownames(seurat_ordered)]

common <- intersect(mes_valid, epi_valid)
mes_valid <- setdiff(mes_valid, common)
epi_valid <- setdiff(epi_valid, common)

gene_list <- list(
  Epithelial_Score = epi_valid,
  Mesenchymal_Score = mes_valid
)

seurat_ordered <- AddModuleScore(
  object = seurat_ordered,
  features = gene_list,
  name = "Enrichment_score" 
)

seurat_ordered$Epithelial_Enrichment <- seurat_ordered$Enrichment_score1
seurat_ordered$Mesenchymal_Enrichment <- seurat_ordered$Enrichment_score2

seurat_ordered$EMT_Difference <- seurat_ordered$Mesenchymal_Enrichment - seurat_ordered$Epithelial_Enrichment

p_diff <- FeaturePlot(seurat_ordered, 
                      reduction = "umap_wnn", 
                      features = "EMT_Difference", 
                      order = TRUE, 
                      pt.size = 0)
print(p_diff)

#remove the axis labels and center is by default (check what is this default value) 
p_diff_final <- p_diff + 
  scale_colour_gradientn(colours = c("blue", "white", "red")) +
  coord_fixed() + 
  theme_minimal() +
  theme(
    legend.position = "right", 
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    aspect.ratio = 1,
    plot.margin = margin(0, 0, 0, 0, "cm"),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(title = "Epithelial - Mesenchymal Score",
       colour = "Difference") +
  scale_x_continuous(expand = c(0.01, 0.01)) + 
  scale_y_continuous(expand = c(0.01, 0.01))

print(p_diff_final)

#remove the axis labels and forcefully set the center/mid point to be at 0
p_diff_final <- p_diff + 
  scale_colour_gradient2(
    low = "blue", 
    mid = "white", 
    high = "red", 
    midpoint = 0  # This is the key line
  ) +
  coord_fixed() + 
  theme_minimal() +
  theme(
    legend.position = "right", 
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    aspect.ratio = 1,
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(title = "Epithelial - Mesenchymal Score", color = "Difference")
print(p_diff_final)


features_to_plot <- c("Epithelial_Enrichment", "Mesenchymal_Enrichment")

plot_df <- seurat_ordered@meta.data %>%
  mutate(sample = Idents(seurat_ordered)) %>%
  select(sample, Epithelial_Enrichment, Mesenchymal_Enrichment)

write.table(plot_df, 
            file = "Epi_mes_listC_enrichment_scores_mut_time_points.txt", 
            sep = "\t", 
            quote = FALSE, 
            row.names = TRUE)

epi_median <- median(seurat_ordered$Epithelial_Enrichment, na.rm = TRUE)
mes_median <- median(seurat_ordered$Mesenchymal_Enrichment, na.rm = TRUE)


clean_umap <- function(p, midpoint_val, title_str) {
  p + scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = midpoint_val) +
    coord_fixed() + 
    theme_minimal() +
    theme(
      legend.position = "right", 
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      aspect.ratio = 1,
      # These four lines remove the extra "padding" around the plot
      plot.margin = margin(0, 0, 0, 0, "cm"),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    ) +
    labs(title = title_str) +
    # This forces the axis to wrap tightly around the data points
    scale_x_continuous(expand = c(0.01, 0.01)) + 
    scale_y_continuous(expand = c(0.01, 0.01))
}

# 3. Generate the plots
p1 <- FeaturePlot(seurat_ordered, reduction = "umap_wnn", features = "Epithelial_Enrichment", order = TRUE, pt.size = 0)
p2 <- FeaturePlot(seurat_ordered, reduction = "umap_wnn", features = "Mesenchymal_Enrichment", order = TRUE, pt.size = 0)

p_epi_final <- clean_umap(p1, epi_median, "Epithelial Score")
p_mes_final <- clean_umap(p2, mes_median, "Mesenchymal Score")

combined_plot <- p_epi_final + p_mes_final
print(combined_plot)
# SAVE THE IMAGE: 
# 300 DPI is standard for publication
# We adjust width/height to fit the square plots perfectly
ggsave("epi_mes_cell_line_score_wnnumap_mut_time_points.png", 
       plot = combined_plot, 
       width = 10,    # Total width for two squares side-by-side
       height = 5,    # Height of one square
       dpi = 1200, 
       bg = "white") # Ensures no transparency issues







