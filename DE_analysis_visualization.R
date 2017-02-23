rm(list = ls())
options <- commandArgs(trailingOnly = TRUE)
print(sprintf("source code is %s", as.character(options[1])))
source(as.character(options[1]))
####

resDir <- as.character(options[2])
print("----")
workDir <- resDir
setwd(workDir)

group1 <- as.character(options[3])
group2 <- as.character(options[4])

samp.cutoff <- as.numeric(options[5])
cpm.val <- as.numeric(options[6])

padj.val <- as.numeric(options[7])
fc.val <- as.numeric(options[8])

####
# source("/Volumes/yli/DE_tools/DE_analysis_vis_Fn.R")
# resDir <- "/Volumes/yli/DE_tools/DEG_analysis_res/"
# workDir <- resDir
# setwd(workDir)
# group1 <- "Control"
# group2 <- "Control_Hypothroy"
# 
# samp.cutoff <- 2
# cpm.val <- 3
# 
# padj.val <- 0.05
# fc.val <- 1.5
####
library(ggplot2)
print("-----")
edgeR.res.fname     <- paste(resDir, "/edgeR-res-comp-", as.character(group1), "_", as.character(group2), "/edgeR-full-GLMDisp-est-cpmCutoof_", cpm.val, "-sampCutoff_", samp.cutoff, ".txt", sep = "")
deseq2.res.fname    <- paste(resDir, "/DEseq2-res-comp-", as.character(group1), "_", as.character(group2), "/DEseq2-full-GLMDisp-est-cpmCutoof_", cpm.val, "-sampCutoff_", samp.cutoff, ".txt", sep = "")
limmavoom.res.fname <- paste(resDir, "/limma-voom-res-comp-", as.character(group1), "_", as.character(group2), "/limma-voom-full-GLMDisp-est-cpmCutoof_", cpm.val, "-sampCutoff_", samp.cutoff, ".txt", sep = "")

plotDir <- paste(as.character(workDir), "/plot-res", sep = "")
if (!dir.exists(plotDir)) dir.create(plotDir)

##Volcano plot
volcano.plot.fn.edger(input = edgeR.res.fname, fname = paste(plotDir, "/edgeR-volcano.pdf", sep = ""), pval = padj.val, fc = fc.val)
volcano.plot.fn.deseq2(input = deseq2.res.fname, fname = paste(plotDir, "/DESeq2-volcano.pdf", sep = ""), pval = padj.val, fc = fc.val)
volcano.plot.fn.voom(input = limmavoom.res.fname, fname = paste(plotDir, "/Voom-volcano.pdf", sep = ""), pval = padj.val, fc = fc.val)

res.edgeR.deg <- get.deg(edgeR.res.fname, fc = fc.val, pval = padj.val, method = 'edgeR')
res.deseq2.deg <- get.deg(deseq2.res.fname, fc = fc.val, pval = padj.val, method = 'DESeq2')
res.voom.deg <- get.deg(limmavoom.res.fname, fc = fc.val, pval = padj.val, method = 'voom')
print(sprintf("No. of DEGs from edgeR = %s, DESeq2 = %s, Limma-voom = %s", 
              length(rownames(res.edgeR.deg)), length(rownames(res.deseq2.deg)), length(rownames(res.voom.deg))))
library(gplots)
plot.name <- paste(plotDir, "/venn.pdf", sep = "")
ol <- venn(list('edgeR' = rownames(res.edgeR.deg),
                'DESeq2' = rownames(res.deseq2.deg),
                'Limma-voom' = rownames(res.voom.deg)), show.plot = F)
pdf(plot.name, height = 5, width = 5)
par(mar = c(1, 0, 0, 0))
plot(ol)
dev.off()
print("-----")
