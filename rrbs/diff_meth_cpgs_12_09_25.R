#the script/code is to get the CpGs; separate code to get the DMRs
#BiocManager::install("methylKit")
library('methylKit')

#sample.id <- list("soft_selected_rep1", "soft_selected_rep2", "soft_selected_rep3", "soft_selected_rep4",
#          "soft_ancestral_rep1", "soft_ancestral_rep2","soft_ancestral_rep3","soft_ancestral_rep4")
sample.id <- list("stiff_selected_rep1", "stiff_selected_rep2", "stiff_selected_rep3", "stiff_selected_rep4",
                  "stiff_ancestral_rep1", "stiff_ancestral_rep2","stiff_ancestral_rep3","stiff_ancestral_rep4")
#sample.id <- list("soft_ancestral_rep1", "soft_ancestral_rep2","soft_ancestral_rep3", "soft_ancestral_rep4",
#                  "stiff_ancestral_rep1", "stiff_ancestral_rep2","stiff_ancestral_rep3", "stiff_ancestral_rep4") #to run for soft ancestral versus stiff ancestral comparison
treatment <- c(1, 1, 1, 1, 0, 0 ,0 ,0)  # 0 = control, 1 = treatment

myobj <- methRead(list("MDA-MB-231-adapted-308Kpa-R4_S9_L001_R1_001_trimmed_bismark_bt2.bismark.cov.gz",
                       "MDA-MB-231-adapted-308Kpa-R4_S9_L001_R1_001_trimmed_bismark_bt2.bismark.cov.gz",
                       "MDA-MB-231-adapted-308Kpa-R2_S10_L001_R1_001_trimmed_bismark_bt2.bismark.cov.gz",
                       "MDA-MB-231-adapted-308Kpa-R1_S7_L001_R1_001_trimmed_bismark_bt2.bismark.cov.gz",
                      "MDA-MB-231-ancestral-308Kpa-R4_S12_L002_R1_001_trimmed_bismark_bt2.bismark.cov.gz",
                       "MDA-MB-231-ancestral-308Kpa-R3_S14_L002_R1_001_trimmed_bismark_bt2.bismark.cov.gz",
                       "MDA-MB-231-ancestral-308Kpa-R2_S13_L002_R1_001_trimmed_bismark_bt2.bismark.cov.gz",
                       "MDA-MB-231-ancestral-308Kpa-R1_S15_L002_R1_001_trimmed_bismark_bt2.bismark.cov.gz"),
                  sample.id=sample.id,
                  pipeline = "bismarkCoverage",
                  assembly="hg38",
                  treatment=treatment,
                  mincov = 10,
                  context = "CpG"
)
head(myobj[[4]])

#get the descriptive stats
getMethylationStats(myobj[[4]], plot=TRUE, both.strands=FALSE)
getCoverageStats(myobj[[4]], plot=TRUE, both.strands=FALSE)

myobj.filt <- filterByCoverage(myobj,
                               lo.count=10,
                               lo.perc=NULL,
                               hi.count=NULL,
                               hi.perc=99.9)
myobj.filt.norm <- normalizeCoverage(myobj.filt, method = "median")
meth <- unite(myobj.filt.norm, destrand=FALSE)

#getCorrelation(meth,plot=TRUE)
#clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)

#PCASamples(meth)

#a filtering step is required to remove sites with little or no variation and SNPs 
#not doing 
myDiff <- calculateDiffMeth(meth,
                            overdispersion = "MN",
                            adjust="BH")

myDiff_df <- getData(myDiff)
write.table(myDiff_df,
            file = "st_sel_st_anc_all_cpgs.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)


########### run the script till here to get all the diff meth cpgs or dmrs ##################
# Simple volcano plot to get an overview of differential methylation
plot(myDiff$meth.diff, -log10(myDiff$qvalue))
abline(v=0)
diffMethPerChr(myDiff)

# get hyper methylated bases and order by qvalue
myDiff25p.hyper <- getMethylDiff(myDiff,
                                 difference=25,
                                 qvalue=0.01,
                                 type="hyper")
myDiff25p.hyper <- myDiff25p.hyper[order(myDiff25p.hyper$qvalue),]

# get hypo methylated bases and order by qvalue
myDiff25p.hypo <- getMethylDiff(myDiff,
                                difference=25,
                                qvalue=0.01,
                                type="hypo")
myDiff25p.hypo <- myDiff25p.hypo[order(myDiff25p.hypo$qvalue),]

# get all differentially methylated bases and order by qvalue
myDiff25p <- getMethylDiff(myDiff,
                           difference=25,
                           qvalue=0.01)
myDiff25p <- myDiff25p[order(myDiff25p$qvalue),]
View(myDiff25p.hyper)

#save the hyper-methylated sites
write.table(myDiff25p.hyper,
            file = "myDiff25p_hyper.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

# Save hypo-methylated sites
write.table(myDiff25p.hypo,
            file = "myDiff25p_hypo.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

cpg_anot <- readFeatureFlank("E:\\tanmay_lele_19thmay\\RRBS\\fastq\\bismark_meth_extractor\\coverage_files\\hg38_cpg_islands_ucsc_12_08.txt", feature.flank.name = c("CpGi", "shores"), flank=2000)
diffCpGann <- annotateWithFeatureFlank(as(myDiff25p,"GRanges"), feature = cpg_anot$CpGi, flank = cpg_anot$shores, feature.name = "CpGi", flank.name = "shores")

# See whether the CpG in myDiff25p belong to a CpG Island or Shore
head(getMembers(diffCpGann))

# This can also be summarized for all differentially methylated CpGs
plotTargetAnnotation(diffCpGann, main = "Differential Methylation Annotation")   



