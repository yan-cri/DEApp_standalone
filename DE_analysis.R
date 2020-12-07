rm(list = ls())
options <- commandArgs(trailingOnly = TRUE)

print(sprintf("source code is %s", paste(getwd(), 'R', 'DE-Analysis-AllComb-Fn.R', sep = '/')))
source(paste(getwd(), 'R', 'DE-Analysis-AllComb-Fn.R', sep = '/'))
de.vis.fn.fname <- paste(getwd(), 'R', 'DE_analysis_vis_Fn.R', sep = '/')

countTabpath <- as.character(options[1])
metaTabpath <- as.character(options[2])
group1 <- as.character(options[3])
group2 <- as.character(options[4])
outputDir <- as.character(options[5])
if (!dir.exists(outputDir)) dir.create(outputDir)
samp.cutoff <- as.numeric(options[6])
cpm.val <- as.numeric(options[7])
saveRes <- as.logical(options[8])
padjOn <- as.logical(options[9])
padj.val <- as.numeric(options[10])
fc.val <- as.numeric(options[11])

# countTabpath = paste(getwd(), 'testData', 'pnas-count_singleFactor.txt', sep = '/')
# metaTabpath = paste(getwd(), 'testData', 'pnas-count_singleFactor-meta.txt', sep = '/')
# cpm.val <- as.numeric(1)
# samp.cutoff <- as.numeric(2)
# outputDir <- '/Users/yanli/Desktop/DEApp_standalone/res-test'
# if (!dir.exists(outputDir)) dir.create(outputDir)
# saveRes <- as.logical(T)
# padjOn <- as.logical(T)
# cpmCutoff = as.numeric(cpm.val)
# sampCutoff=as.numeric(samp.cutoff)
# group1="DHT"
# group2="Control"
# padjust=as.character(padjOn)
# padj.val <- 0.05
# fc.val <- 1.5

#######################################
###edgeR analysis 
setwd(outputDir)
print(sprintf('currect work directory is %s', as.character(outputDir)))
print("============================")
print("Start edgeR DE analysis")
edgeR.res.path <- paste(getwd(), "/edgeR-res-comp-", as.character(trim(group1)), "_", as.character(trim(group2)), sep="")
##This above path need to be kept the same in the source code, "/group/bioinformatics/Pipelines/Development/RNAseq_beta/RNAseq.v9.0.1/build/module-ve/DE-Analysis-AllComb.R"
print(sprintf("Raw count is from %s, low expression is removed at cpm = %.0f, group no = % .0f", 
              as.character(countTabpath),
              as.numeric(cpm.val),
              as.numeric(samp.cutoff)))
GeneRegFilter="both"
print(sprintf("The filter is on %s", as.character(as.character(GeneRegFilter))))
print(sprintf("The comparision is between %s and %s", as.character(trim(group1)), as.character(trim(group2))))
##The FC is calculated on group2/group1
print(sprintf("The results are saved in %s", as.character(edgeR.res.path)))
edgeR.res <- edgeR.de.analysis(countTabpath = as.character(countTabpath), 
                               metaTabpath = as.character(metaTabpath), 
                               cpmCutoff = as.numeric(cpm.val), sampCutoff=as.numeric(samp.cutoff), 
                               group1=as.character(trim(group1)), group2=as.character(trim(group2)), 
                               padjust=as.character(padjOn), Pcutoff=as.numeric(padj.val), FCcutoff=as.numeric(fc.val), 
                               GeneRegFilter=as.character(GeneRegFilter), 
                               saveRes = as.character(saveRes))

if(padjOn){
  print(sprintf("edgeR identified %s DEGs at cpm = %s, samp = %s (total features %s), FDR adj-p = %s, and FC = %s.", 
                length(edgeR.res$degName),
                as.numeric(cpm.val), 
                as.numeric(samp.cutoff), 
                dim(edgeR.res$res)[1], 
                as.numeric(padj.val),
                as.numeric(fc.val)) )
}else{
  print(sprintf("edgeR identified %s DEGs at cpm = %s, samp = %s (total features %s), original p-value = %s, and FC = %s.", 
                length(edgeR.res$degName),
                as.numeric(cpm.val), 
                as.numeric(samp.cutoff), 
                dim(edgeR.res$res)[1], 
                as.numeric(padj.val),
                as.numeric(fc.val)) )
}

print("**********")
print("End edgeR DE analysis")
print("============================")

print("============================")
print("Start DESeq2 DE analysis")
deseq2.res.path <- paste(getwd(), "/DEseq2-res-comp-", as.character(trim(trim(group1))), "_", as.character(trim(group2)), sep="")
print(sprintf("The results are saved in %s", as.character(deseq2.res.path)))
deseq2.res <- deseq2.de.analysis(countTabpath = as.character(countTabpath),
                                 metaTabpath = as.character(metaTabpath), 
                                 cpmCutoff = as.numeric(cpm.val), sampCutoff=as.numeric(samp.cutoff), 
                                 group1=as.character(trim(group1)), group2=as.character(trim(group2)),
                                 padjust=as.character(padjOn), Pcutoff=as.numeric(padj.val), FCcutoff=as.numeric(fc.val), 
                                 GeneRegFilter=as.character(GeneRegFilter), 
                                 saveRes = as.character(saveRes))
if(padjOn){
  print(sprintf("deseq2 identified %s DEGs at cpm = %s, samp = %s (total features %s), FDR adj-p = %s, and FC = %s.", 
                length(deseq2.res$degName),
                as.numeric(cpm.val), 
                as.numeric(samp.cutoff), 
                dim(deseq2.res$res)[1], 
                as.numeric(padj.val),
                as.numeric(fc.val)) )
}else{
  print(sprintf("deseq2 identified %s DEGs at cpm = %s, samp = %s (total features %s), original p-value = %s, and FC = %s.", 
                length(deseq2.res$degName),
                as.numeric(cpm.val), 
                as.numeric(samp.cutoff), 
                dim(deseq2.res$res)[1], 
                as.numeric(padj.val),
                as.numeric(fc.val)) )
}

print("**********")
print("End DESeq2 DE analysis")
print("============================")
print("============================")
print("Start Limma-Voom DE analysis")
limmavoom.res.path <- paste(getwd(), "/limma-voom-res-comp-", as.character(trim(group1)), "_", as.character(trim(group2)), sep="")
##This above path need to be kept the same in the source code, "/group/bioinformatics/Pipelines/Development/RNAseq_beta/RNAseq.v9.0.1/build/module-ve/DE-Analysis-AllComb.R"
print(sprintf("The results are saved in %s", as.character(limmavoom.res.path)))
limma.res <- limmavoom.de.analysis(countTabpath = as.character(countTabpath), 
                                   metaTabpath = as.character(metaTabpath), 
                                   cpmCutoff = as.numeric(cpm.val), sampCutoff=as.numeric(samp.cutoff), 
                                   group1=as.character(trim(group1)), group2=as.character(trim(group2)), 
                                   padjust=as.character(padjOn), Pcutoff=as.numeric(padj.val), FCcutoff=as.numeric(fc.val), 
                                   GeneRegFilter=as.character(GeneRegFilter), 
                                   saveRes = as.character(saveRes))
if(padjOn){
  print(sprintf("limma identified %s DEGs at cpm = %s, samp = %s (total features %s), FDR adj-p = %s, and FC = %s.", 
                length(limma.res$degName),
                as.numeric(cpm.val), 
                as.numeric(samp.cutoff), 
                dim(limma.res$res)[1], 
                as.numeric(padj.val),
                as.numeric(fc.val)) )
}else{
  print(sprintf("limma identified %s DEGs at cpm = %s, samp = %s (total features %s), original p-value = %s, and FC = %s.", 
                length(limma.res$degName),
                as.numeric(cpm.val), 
                as.numeric(samp.cutoff), 
                dim(limma.res$res)[1], 
                as.numeric(padj.val),
                as.numeric(fc.val)) )
}

print("**********")
print("End Limma-Voom DE analysis")
print("============================")
########################################
########################################
########################################
print("============================")
print("Start DE analysis visualization analysis")
source(as.character(de.vis.fn.fname))
print(sprintf("DE visualization function source code is %s", as.character(de.vis.fn.fname) ))

print("-----")
print('Start volcano plots')
library(ggplot2)
print("-----")
edgeR.res.fname     <- paste(outputDir, "/edgeR-res-comp-", as.character(group1), "_", as.character(group2), "/edgeR-full-GLMDisp-est-cpmCutoof_", cpm.val, "-sampCutoff_", samp.cutoff, ".txt", sep = "")
deseq2.res.fname    <- paste(outputDir, "/DEseq2-res-comp-", as.character(group1), "_", as.character(group2), "/DEseq2-full-GLMDisp-est-cpmCutoof_", cpm.val, "-sampCutoff_", samp.cutoff, ".txt", sep = "")
limmavoom.res.fname <- paste(outputDir, "/limma-voom-res-comp-", as.character(group1), "_", as.character(group2), "/limma-voom-full-GLMDisp-est-cpmCutoof_", cpm.val, "-sampCutoff_", samp.cutoff, ".txt", sep = "")

plotDir <- paste(as.character(outputDir), "/plot-res", sep = "")
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
print('End volcano plots')
print("-----")
###Venn part
print("-----")
print('Start DEGs overlapping analysis')
library(gplots)
olDir <- paste(as.character(outputDir), "/ol-res", sep = "")
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

print('End DEGs overlapping analysis')
print("-----")
print("End DE analysis visualization analysis")
print("============================")

