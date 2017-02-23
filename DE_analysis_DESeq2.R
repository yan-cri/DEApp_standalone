rm(list = ls())
source("/group/bioinformatics/Pipelines/Development/RNAseq_beta/RNAseq.v9.0.1/build/module-ve/DE-Analysis-AllComb.R")
options <- commandArgs(trailingOnly = TRUE)
cpm.val <- as.numeric(options[7])
samp.cutoff <- as.numeric(options[6])
padj.val <- c(0.01, 0.05, 0.1)
fc.val <- c(0, 1.5, 2)
saveRes <- as.logical(options[8])
padjOn <- as.logical(options[9])

#######################################
###########################################
###DEseq2 analysis 
setwd(as.character(options[5]))
print("============================")
print("Start DESeq2 DE analysis")
deseq2.res.path <- paste(getwd(), "/DEseq2-res-comp-", as.character(trim(options[3])), "_", as.character(trim(options[4])), sep="")
##This above path need to be kept the same in the source code, "/group/bioinformatics/Pipelines/Development/RNAseq_beta/RNAseq.v9.0.1/build/module-ve/DE-Analysis-AllComb.R"
print(sprintf("Raw count is from %s, low expression is removed at cpm = %.0f, group no = % .0f", 
              as.character(options[1]),
              as.numeric(cpm.val),
              as.numeric(samp.cutoff)))
GeneRegFilter="both"
print(sprintf("The filter is on %s", as.character(as.character(GeneRegFilter))))
print(sprintf("The comparision is between %s and %s", as.character(trim(options[3])), as.character(trim(options[4]))))
##The FC is calculated on group2/group1
print(sprintf("The results are saved in %s", as.character(deseq2.res.path)))
de.num.padj.fc.deseq2 <- matrix(NA, nrow = 3, ncol = 3)
for (i in 1:3){
  for (j in 1:3){
    deseq2.res <- deseq2.de.analysis(countTabpath = as.character(options[1]),
                                     metaTabpath = as.character(options[2]), 
                                     cpmCutoff = as.numeric(cpm.val), sampCutoff=as.numeric(samp.cutoff), 
                                     group1=as.character(options[3]), group2=as.character(options[4]),
                                     padjust=as.character(padjOn), Pcutoff=as.numeric(padj.val[j]), FCcutoff=as.numeric(fc.val[i]), 
                                     GeneRegFilter=as.character(GeneRegFilter), 
                                     saveRes = as.character(saveRes))
    print(sprintf("deseq2 identified %s DEGs at cpm = %s, samp = %s (total features %s), FDR adj-p = %s, and FC = %s.", 
                  length(deseq2.res$degName),
                  as.numeric(cpm.val), 
                  as.numeric(samp.cutoff), 
                  dim(deseq2.res$res)[1], 
                  as.numeric(padj.val[j]),
                  as.numeric(fc.val[i])) )
    de.num.padj.fc.deseq2[i,j] = length(deseq2.res$degName)
  }
}
colnames(de.num.padj.fc.deseq2) <- paste("adj-p < ", as.numeric(padj.val), sep = "")
rownames(de.num.padj.fc.deseq2) <- paste("|FC| > ", as.numeric(fc.val), sep = "")
fname <- paste(as.character(deseq2.res.path),"/DESeq2-summary-filter-", as.character(GeneRegFilter), ".txt", sep = "")
write.table(de.num.padj.fc.deseq2, as.character(fname), quote = F, col.names = NA, row.names = T, sep = "\t")
print("**********")
GeneRegFilter="up"
print(sprintf("The filter is on %s", as.character(as.character(GeneRegFilter))))
de.num.padj.fc.deseq2 <- matrix(NA, nrow = 3, ncol = 3)
for (i in 1:3){
  for (j in 1:3){
    deseq2.res <- deseq2.de.analysis(countTabpath = as.character(options[1]),
                                     metaTabpath = as.character(options[2]), 
                                     cpmCutoff = as.numeric(cpm.val), sampCutoff=as.numeric(samp.cutoff), 
                                     group1=as.character(options[3]), group2=as.character(options[4]),
                                     padjust=as.character(padjOn), Pcutoff=as.numeric(padj.val[j]), FCcutoff=as.numeric(fc.val[i]), 
                                     GeneRegFilter=as.character(GeneRegFilter), 
                                     saveRes = as.character(saveRes))
    print(sprintf("deseq2 identified %s DEGs at cpm = %s, samp = %s (total features %s), FDR adj-p = %s, and FC = %s.", 
                  length(deseq2.res$degName),
                  as.numeric(cpm.val), 
                  as.numeric(samp.cutoff), 
                  dim(deseq2.res$res)[1], 
                  as.numeric(padj.val[j]),
                  as.numeric(fc.val[i])) )
    de.num.padj.fc.deseq2[i,j] = length(deseq2.res$degName)
  }
}
colnames(de.num.padj.fc.deseq2) <- paste("adj-p < ", as.numeric(padj.val), sep = "")
rownames(de.num.padj.fc.deseq2) <- paste("|FC| > ", as.numeric(fc.val), sep = "")
fname <- paste(as.character(deseq2.res.path),"/DESeq2-summary-filter-", as.character(GeneRegFilter), ".txt", sep = "")
write.table(de.num.padj.fc.deseq2, as.character(fname), quote = F, col.names = NA, row.names = T, sep = "\t")
print("**********")
GeneRegFilter="down"
print(sprintf("The filter is on %s", as.character(as.character(GeneRegFilter))))
de.num.padj.fc.deseq2 <- matrix(NA, nrow = 3, ncol = 3)
for (i in 1:3){
  for (j in 1:3){
    deseq2.res <- deseq2.de.analysis(countTabpath = as.character(options[1]),
                                     metaTabpath = as.character(options[2]), 
                                     cpmCutoff = as.numeric(cpm.val), sampCutoff=as.numeric(samp.cutoff), 
                                     group1=as.character(options[3]), group2=as.character(options[4]),
                                     padjust=as.character(padjOn), Pcutoff=as.numeric(padj.val[j]), FCcutoff=as.numeric(fc.val[i]), 
                                     GeneRegFilter=as.character(GeneRegFilter), 
                                     saveRes = as.character(saveRes))
    print(sprintf("deseq2 identified %s DEGs at cpm = %s, samp = %s (total features %s), FDR adj-p = %s, and FC = %s.", 
                  length(deseq2.res$degName),
                  as.numeric(cpm.val), 
                  as.numeric(samp.cutoff), 
                  dim(deseq2.res$res)[1], 
                  as.numeric(padj.val[j]),
                  as.numeric(fc.val[i])) )
    de.num.padj.fc.deseq2[i,j] = length(deseq2.res$degName)
  }
}
colnames(de.num.padj.fc.deseq2) <- paste("adj-p < ", as.numeric(padj.val), sep = "")
rownames(de.num.padj.fc.deseq2) <- paste("|FC| > ", as.numeric(fc.val), sep = "")
fname <- paste(as.character(deseq2.res.path),"/DESeq2-summary-filter-", as.character(GeneRegFilter), ".txt", sep = "")
write.table(de.num.padj.fc.deseq2, as.character(fname), quote = F, col.names = NA, row.names = T, sep = "\t")
print("End DESeq2 DE analysis")
print("============================")
