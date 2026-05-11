library('methylKit')
library("genomation")
library("GenomicRanges")

sample.id <- list("soft_selected_rep1", "soft_selected_rep2", "soft_selected_rep3", "soft_selected_rep4",
                  "soft_ancestral_rep1", "soft_ancestral_rep2","soft_ancestral_rep3","soft_ancestral_rep4")
#sample.id <- list("soft_ancestral_rep1", "soft_ancestral_rep2","soft_ancestral_rep3", "soft_ancestral_rep4",
 #                 "stiff_ancestral_rep1", "stiff_ancestral_rep2","stiff_ancestral_rep3", "stiff_ancestral_rep4") #to run for soft ancestral versus stiff ancestral comparison
treatment <- c(1,1,1,1,0,0,0,0)  # 0 = control, 1 = treatment
myobj <- methRead(list("MDA-MB-231-adapted-1kPa-R1_S2_L001_R1_001_trimmed_bismark_bt2.bismark.cov.gz",
                       "MDA-MB-231-adapted-1kPa-R2_S5_L001_R1_001_trimmed_bismark_bt2.bismark.cov.gz",
                       "MDA-MB-231-adapted-1kPa-R3_S11_L002_R1_001_trimmed_bismark_bt2.bismark.cov.gz",
                       "MDA-MB-231-adapted-1kPa-R4_S4_L001_R1_001_trimmed_bismark_bt2.bismark.cov.gz",
                       "MDA-MB-231-ancestral-1kPa-R1_S17_L002_R1_001_trimmed_bismark_bt2.bismark.cov.gz",
                       "MDA-MB-231-ancestral-1kPa-R2_S18_L002_R1_001_trimmed_bismark_bt2.bismark.cov.gz",
                       "MDA-MB-231-ancestral-1kPa-R3_S19_L002_R1_001_trimmed_bismark_bt2.bismark.cov.gz",
                       "MDA-MB-231-ancestral-1kPa-R4_S16_L002_R1_001_trimmed_bismark_bt2.bismark.cov.gz"),
                  sample.id=sample.id,
                  pipeline = "bismarkCoverage",
                  assembly="hg38",
                  treatment=treatment,
                  mincov = 10,
                  context = "CpG"
)
myobj.filt <- filterByCoverage(myobj,
                               lo.count=10,
                               lo.perc=NULL,
                               hi.count=NULL,
                               hi.perc=99.9)
myobj.filt.norm <- normalizeCoverage(myobj.filt, method = "median")

tiles <- tileMethylCounts(myobj.filt.norm,              # from your previous filtering + normalization
                          win.size=1000,                # 1000 bp windows
                          step.size=1000,               # non-overlapping
                          cov.bases=3,                  # min CpG bases per window
                          mc.cores=4)

tiles_united <- unite(tiles, destrand=FALSE)
dmr_tiles <- calculateDiffMeth(tiles_united, overdispersion="MN", adjust="BH")

#save all these dmr tiles
dmr_tiles_df <- getData(dmr_tiles)

# Save to tab-separated text file
write.table(dmr_tiles_df,
            file="so_sel_so_anc_all_dmrs_allreps.txt",
            sep="\t",
            row.names=FALSE,
            quote=FALSE)







dmrs_25p <- getMethylDiff(dmr_tiles,
                          difference = 25,
                          qvalue = 0.05)
dmrs_25p <- dmrs_25p[order(dmrs_25p$qvalue), ]

# Hyper-methylated regions
dmrs_25p_hyper <- getMethylDiff(dmr_tiles,
                                difference = 25,
                                qvalue = 0.05,
                                type = "hyper")
dmrs_25p_hyper <- dmrs_25p_hyper[order(dmrs_25p_hyper$qvalue), ]

# Hypo-methylated regions
dmrs_25p_hypo <- getMethylDiff(dmr_tiles,
                               difference = 25,
                               qvalue = 0.05,
                               type = "hypo")
dmrs_25p_hypo <- dmrs_25p_hypo[order(dmrs_25p_hypo$qvalue), ]

# 💾 Step 5: Save DMRs as tab-separated text files
write.table(dmrs_25p, "DMRs_all_q05_diff25.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(dmrs_25p_hyper, "DMRs_hyper_q05_diff25.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(dmrs_25p_hypo, "DMRs_hypo_q05_diff25.txt", sep = "\t", quote = FALSE, row.names = FALSE)