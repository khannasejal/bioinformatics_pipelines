library(DiffBind)

#read the sample sheet metadata
samples <- read.delim("samples.txt", sep = "\t", header = TRUE)
View(samples)
dbObj <- dba(sampleSheet = samples)

dbObj <- dba.count(dbObj)
dbObj <- dba.contrast(dbObj, categories=DBA_CONDITION)
dbObj <- dba.analyze(dbObj)
report <- dba.report(dbObj)
write.csv(as.data.frame(report), "DiffBind_results_selected_samples_b1.csv")

dba.plotPCA(
  dbObj,
  attributes = DBA_CONDITION,   # color by your Condition column
  label = DBA_ID,               # label using SampleID
  bLog = TRUE                   # log-transform (recommended for ATAC-seq)
)