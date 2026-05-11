library(COCONUT)
base_dir_intensities <- "D:/sejal/RNA_seq_dec24th/series_matrix/GSE_folder_cocunut/dec30th_more_datasets/paste_folder/RNA/new_march15th/intersect"
base_dir_metadata <- "D:/sejal/RNA_seq_dec24th/series_matrix/GSE_folder_cocunut/dec30th_more_datasets/paste_folder_metadata/RNA"

gse_files <- list.files(base_dir_intensities, pattern = "_intensities_coconut_intersect.txt", full.names = TRUE)
#gse_ids <- unique(sub("_.*", "", basename(gse_files)))
gse_ids <- c(
  "GSE100075","GSE118713", "GSE103243", "GSE104872", "GSE104985", "GSE106681", "GSE108787",
  "GSE115609","GSE117942", "GSE120756", "GSE123285", "GSE126003", "GSE128458", "GSE137270", "GSE140186", "GSE140758", "GSE144378", "GSE147745",
  "GSE149428", "GSE149858", "GSE150997", "GSE152360", "GSE165571", "GSE178303",
  "GSE181909", "GSE182631", "GSE186679", "GSE186682", "GSE190384", "GSE206724",
  "GSE222367", "GSE234074", "GSE241764", "GSE243150", "GSE243454", "GSE247138",
  "GSE251644","GSE253717","GSE266932","GSE78199")
  
#didn't run for:: ("GSE115737","GSE59536","GSE62613") only 1 replicate for control

GSEs <- list()
for (gse_id in gse_ids) {
  print(gse_id)
  signal_file <- file.path(base_dir_intensities, paste0(gse_id, "_intensities_coconut_intersect.txt"))
  print(signal_file)
  metadata_file <- file.path(base_dir_metadata, paste0(gse_id, "_SraRunTable.csv"))
  print(metadata_file)
  
  genes <- as.matrix(read.table(signal_file, sep = "\t", header = TRUE, row.names = 1))
  pheno <- read.csv(metadata_file, sep = ",", header = TRUE, row.names = 1)
  
  GSEs[[gse_id]] <- list(
    pheno = pheno,
    genes = genes   
  )
}

for (gse in GSEs) {
  print(colnames(gse$genes))
  print(rownames(gse$pheno))
}

GSEs_COCONUT <- COCONUT(
  GSEs = GSEs,                         # List of datasets
  control.0.col = "Condition",         # Use the correct column name for control vs treatment
  byPlatform = FALSE                   # Whether to normalize by platform (set to FALSE)
)
gse_ids <- names(GSEs_COCONUT$controlList$GSEs)

for (gse_id in gse_ids) {
  # Combine control and normalized genes for the current GSE ID
  combined_genes <- cbind(
    GSEs_COCONUT$controlList$GSEs[[gse_id]]$genes,   # Control genes
    GSEs_COCONUT$COCONUTList[[gse_id]]$genes         # Normalized genes
  )
  
  # Create a file name
  file_name <- paste0(gse_id, "_coconut_normalised_combined_genes_march15th.txt")
  
  # Save as a text file
  write.table(
    combined_genes,
    file = file_name,
    sep = "\t",
    col.names = NA,    # Include column names
    quote = FALSE      # Avoid quotes around data
  )
  
  # Print a message to indicate progress
  cat("Saved:", file_name, "\n")
}
