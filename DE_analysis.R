rm(list = ls())
options <- commandArgs(trailingOnly = TRUE)

print(sprintf("source code is %s", paste(getwd(), 'R', 'DE-Analysis-AllComb-Fn.R', sep = '/')))
source(paste(getwd(), 'R', 'DE-Analysis-AllComb-Fn.R', sep = '/'))

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
# outputDir <- '/Users/yanli/Desktop/DE_tools/res-test'
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
