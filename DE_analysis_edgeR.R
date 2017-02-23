rm(list = ls())
options <- commandArgs(trailingOnly = TRUE)
print(sprintf("source code is %s", as.character(options[10])))
source(as.character(options[10]))
countTabpath <- as.character(options[1])
metaTabpath <- as.character(options[2])
group1 <- as.character(options[3])
group2 <- as.character(options[4])
workDir <- as.character(options[5])
samp.cutoff <- as.numeric(options[6])
cpm.val <- as.numeric(options[7])
saveRes <- as.logical(options[8])
padjOn <- as.logical(options[9])
padj.val <- as.numeric(options[11])
fc.val <- as.numeric(options[12])

#######################################
########################################
### Based on STAR alignment results
###########################################
###edgeR analysis 
setwd(workDir)
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
for (i in 1:length(fc.val)){
  for (j in 1:length(padj.val)){
    edgeR.res <- edgeR.de.analysis(countTabpath = as.character(countTabpath), 
                                   metaTabpath = as.character(metaTabpath), 
                                   cpmCutoff = as.numeric(cpm.val), sampCutoff=as.numeric(samp.cutoff), 
                                   group1=as.character(trim(group1)), group2=as.character(trim(group2)), 
                                   padjust=as.character(padjOn), Pcutoff=as.numeric(padj.val[j]), FCcutoff=as.numeric(fc.val[i]), 
                                   GeneRegFilter=as.character(GeneRegFilter), 
                                   saveRes = as.character(saveRes))
    print(sprintf("edgeR identified %s DEGs at cpm = %s, samp = %s (total features %s), FDR adj-p = %s, and FC = %s.", 
                  length(edgeR.res$degName),
                  as.numeric(cpm.val), 
                  as.numeric(samp.cutoff), 
                  dim(edgeR.res$res)[1], 
                  as.numeric(padj.val[j]),
                  as.numeric(fc.val[i])) )
    de.num.padj.fc.edger[i,j] = length(edgeR.res$degName)
  }
}
print("End edgeR DE analysis")
print("============================")
