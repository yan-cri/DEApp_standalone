rm(list = ls())
options <- commandArgs(trailingOnly = TRUE)
print(sprintf("source code is %s", paste(getwd(), 'R', 'DE_analysis_vis_Fn.R', sep = '/')))
source(paste(getwd(), 'R', 'DE_analysis_vis_Fn.R', sep = '/'))
####

resDir <- as.character(options[1])
if (!dir.exists(resDir)) dir.create(resDir)
print("----")
workDir <- resDir
setwd(workDir)
print(sprintf("Working directory is at %s, the same as results saving directory", resDir))

group1 <- as.character(options[2])
group2 <- as.character(options[3])

samp.cutoff <- as.numeric(options[4])
cpm.val <- as.numeric(options[5])

padj.val <- as.numeric(options[6])
fc.val <- as.numeric(options[7])

####
# resDir <- "/Volumes/yli/DE_tools/DEG_analysis_res/"
# workDir <- resDir
# setwd(workDir)
# group1 <- "DBDB_Control"
# group2 <- "DBDB"
# 
# samp.cutoff <- 2
# cpm.val <- 3
# 
# padj.val <- 0.05
# fc.val <- 1
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
print(sprintf("Volcano plots are saved at %s", plotDir))
###Venn part
library(gplots)
olDir <- paste(as.character(workDir), "/ol-res", sep = "")
if (!dir.exists(olDir)) dir.create(olDir)
print(sprintf("overlap results are saved at %s", olDir))
plot.name <- paste(olDir, "/venn.pdf", sep = "")
ol <- venn(list('edgeR' = rownames(res.edgeR.deg),
                'DESeq2' = rownames(res.deseq2.deg),
                'Limma-voom' = rownames(res.voom.deg)), show.plot = F)
pdf(plot.name, height = 5, width = 5)
par(mar = c(1, 0, 0, 0))
plot(ol)
dev.off()

##Write ol results
ol1 <- cbind(res.edgeR.deg[match(attr(ol, "intersections")$`edgeR:DESeq2:Limma-voom`, rownames(res.edgeR.deg)),],
             res.voom.deg[match(attr(ol, "intersections")$`edgeR:DESeq2:Limma-voom`, rownames(res.voom.deg)),],
             res.deseq2.deg[match(attr(ol, "intersections")$`edgeR:DESeq2:Limma-voom`, rownames(res.deseq2.deg)),])
ol2 <- cbind(res.edgeR.deg[match(attr(ol, "intersections")$`edgeR:Limma-voom`, rownames(res.edgeR.deg)),],
             res.voom.deg[match(attr(ol, "intersections")$`edgeR:Limma-voom`, rownames(res.voom.deg)),])
ol3 <- cbind(res.edgeR.deg[match(attr(ol, "intersections")$`edgeR:DESeq2`, rownames(res.edgeR.deg)),],
             res.deseq2.deg[match(attr(ol, "intersections")$`edgeR:DESeq2`, rownames(res.deseq2.deg)),])
ol4 <- cbind(res.voom.deg[match(attr(ol, "intersections")$`DESeq2:Limma-voom`, rownames(res.voom.deg)),],
             res.deseq2.deg[match(attr(ol, "intersections")$`DESeq2:Limma-voom`, rownames(res.deseq2.deg)),])
ol5 <- res.deseq2.deg[match(attr(ol, "intersections")$`DESeq2`, rownames(res.deseq2.deg)),]
ol6 <- res.edgeR.deg[match(attr(ol, "intersections")$`edgeR`, rownames(res.edgeR.deg)),]
ol7 <- res.voom.deg[match(attr(ol, "intersections")$`Limma-voom`, rownames(res.voom.deg)),]

colnames(ol1) <- c(paste("edgeR", colnames(res.edgeR.deg), sep = "_"), paste("voom", colnames(res.voom.deg), sep = "_"), paste("DESeq2", colnames(res.deseq2.deg), sep = "_") )
colnames(ol2) <- c(paste("edgeR", colnames(res.edgeR.deg), sep = "_"), paste("voom", colnames(res.voom.deg), sep = "_"))
colnames(ol3) <- c(paste("edgeR", colnames(res.edgeR.deg), sep = "_"), paste("DESeq2", colnames(res.deseq2.deg), sep = "_"))
colnames(ol4) <- c(paste("voom", colnames(res.voom.deg), sep = "_"), paste("DESeq2", colnames(res.deseq2.deg), sep = "_"))
colnames(ol5) <- colnames(res.deseq2.deg)
colnames(ol6) <- colnames(res.edgeR.deg)
colnames(ol7) <- colnames(res.voom.deg)

ol1.save <- ol1[order(ol1$`edgeR_PValue`), c(1,4,5,7,10,11,15,18,19)]
ol2.save <- ol2[order(ol2$`edgeR_PValue`), c(1,4,5,7,10,11)]
ol3.save <- ol3[order(ol3$`edgeR_PValue`), c(1,4,5,8,11,12)]
ol4.save <- ol4[order(ol4$`voom_P.Value`), c(1,4,5,9,12,13)]
ol5.save <- ol5[order(ol5$`pvalue`), c(2,5,6)]
ol6.save <- ol6[order(ol6$`PValue`), c(1,4,5)]
ol7.save <- ol7[order(ol7$`P.Value`), c(1,4,5)]

write.table(ol1.save, file = paste(olDir, "all3_overlap.txt", sep = "/"), sep ="\t", quote = F, row.names = T, col.names = NA)
write.table(ol2.save, file = paste(olDir, "edger_voom_overlap_only.txt", sep = "/"), sep ="\t", quote = F, row.names = T, col.names = NA)
write.table(ol3.save, file = paste(olDir, "edger_deseq2_overlap_only.txt", sep = "/"), sep ="\t", quote = F, row.names = T, col.names = NA)
write.table(ol4.save, file = paste(olDir, "voom_deseq2_overlap_only.txt", sep = "/"), sep ="\t", quote = F, row.names = T, col.names = NA)
write.table(ol5.save, file = paste(olDir, "deseq2_only.txt", sep = "/"), sep ="\t", quote = F, row.names = T, col.names = NA)
write.table(ol6.save, file = paste(olDir, "edgeR_only.txt", sep = "/"), sep ="\t", quote = F, row.names = T, col.names = NA)
write.table(ol7.save, file = paste(olDir, "voom_only.txt", sep = "/"), sep ="\t", quote = F, row.names = T, col.names = NA)

print("-----")
