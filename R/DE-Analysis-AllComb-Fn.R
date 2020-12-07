library(edgeR)
library(limma)
library(DESeq2)

trim <- function (x) gsub("^\\s+|\\s+$", "", x)

edgeR.de.analysis <- function(countTabpath, metaTabpath, cpmCutoff, sampCutoff, group1, group2, padjust, Pcutoff, FCcutoff, GeneRegFilter, saveRes) {
  count.all <- read.delim(as.character(countTabpath), row.names=1, check.names=FALSE)
  metaTab <- read.delim(as.character(metaTabpath), header = T)
  colnames(metaTab) <- tolower(colnames(metaTab))
  metaTab <- metaTab[match(colnames(count.all), metaTab$sample),]
  
  if (dim(metaTab)[2]>2) {
    groupinfo <- metaTab[,2]
    for (i in 3:length(metaTab[1,])) {
      groupinfo <- paste(groupinfo, metaTab[,i], sep=".")
    }
    Group <- factor(groupinfo)
  } else {
    Group <- factor(metaTab[,2])
  }
  
  dge.count <- DGEList(counts=count.all, group=Group)
  dge.count <- calcNormFactors(dge.count)
  
  ### Remove low expression reads from NA/0 values
  cpm.count <- cpm(dge.count$counts) #cpm normalization
  threshold.value <- as.numeric(sampCutoff)
  keep = rowSums(cpm.count > as.numeric(cpmCutoff)) >=threshold.value  
  
  ##(dge.count.rmlow): DGEList of count rm low expression values.
  dge.count.rmlow <- dge.count[keep,]
  dim(dge.count.rmlow$counts)
  dge.count.rmlow$samples$lib.size <- colSums(dge.count.rmlow$counts)
  
  dge.count.rmlow <- calcNormFactors(dge.count.rmlow,method="TMM")
  log.cpm.rmlow <- cpm(dge.count.rmlow, normalized.lib.sizes=T, log=T)
  
  if (saveRes==TRUE) {
    edgeR.res.path <- paste(getwd(), "/edgeR-res-comp-", as.character(trim(group1)), "_", as.character(trim(group2)), sep="")
    if (!dir.exists(edgeR.res.path)) dir.create(edgeR.res.path)
    logcpm.fname <- paste(as.character(getwd()), "/rmlow-logcpm-", "cpmCutoof_", as.character(cpmCutoff), "-sampCutoff_", as.character(sampCutoff), "-", as.character(dim(log.cpm.rmlow)[1]), ".txt", sep=""  )
    if (!file.exists(logcpm.fname)) write.table(log.cpm.rmlow, as.character(logcpm.fname), row.names=T, col.names=NA, quote=F, sep="\t")
  }
  
  pvalue=as.numeric(Pcutoff)
  FC=as.numeric(FCcutoff)
  gfilter <- as.character(GeneRegFilter)
  
  f <- factor(Group, levels=levels(Group))
  design <- model.matrix(~0+f)
  rownames(design) <- rownames(dge.count.rmlow$samples)
  colnames(design) <- levels(Group)
  
  dge.count.rmlow <- estimateDisp(dge.count.rmlow, design)
  dge.count.rmlow$common.dispersion
  summary(dge.count.rmlow$tagwise.dispersion)
  
  glm.fit <- glmFit(dge.count.rmlow, design)
  
  comp <- makeContrasts(contrasts=paste(as.character(trim(group2)), as.character(trim(group1)), sep="-"), levels=design )
  test <- glmLRT(glm.fit, contrast=comp)
  
  tpsave <- topTags(test,n=Inf, adjust.method="BH", sort.by="PValue")
  
  if (saveRes==TRUE) {
    tp.fname <- paste(as.character(edgeR.res.path), "/edgeR-full-GLMDisp-est-", "cpmCutoof_", as.character(cpmCutoff), "-sampCutoff_", as.character(sampCutoff), ".txt", sep="" )
    # print(head(tpsave$table))
    if (!file.exists(tp.fname)) write.table(tpsave$table, as.character(tp.fname), row.names=T, col.names=NA, quote=F, sep="\t")
  }
  
  if(as.character(padjust)==TRUE) {
    tp <- topTags(test,n=Inf, adjust.method="BH", sort.by="none")
    filter <- decideTestsDGE(test, adjust.method="BH", p.value=pvalue, lfc=log2(FC))
    tp$table$filter <- as.numeric(filter)
  } else {
    tp <- topTags(test,n=Inf, adjust.method="none", sort.by="none")
    filter <- decideTestsDGE(test, adjust.method="none", p.value=pvalue, lfc=log2(FC))
    tp$table$filter <- as.numeric(filter)
  }
  # print(summary(filter))
  
  if (gfilter == "up") {
    upreg <- subset(tp$table, filter==1)
    if (saveRes==TRUE) {
      filter.name <- paste(as.character(edgeR.res.path), "/edgeR-upreg-", "cpmCutoff_", as.character(cpmCutoff), "-sampCutoff_", as.character(sampCutoff), "-fdr_", as.character(pvalue), "-FC_",as.character(FC), "-GeneName.txt", sep="" )
      filter.name.full <- paste(as.character(edgeR.res.path), "/edgeR-upreg-", "cpmCutoff_", as.character(cpmCutoff), "-sampCutoff_", as.character(sampCutoff), "-fdr_", as.character(pvalue), "-FC_",as.character(FC), "-full.txt", sep="" )
    }
  } else if (gfilter == "down") {
    upreg <- subset(tp$table, filter==-1)
    if (saveRes==TRUE){
      filter.name <- paste(as.character(edgeR.res.path), "/edgeR-downreg-", "cpmCutoff_", as.character(cpmCutoff), "-sampCutoff_", as.character(sampCutoff), "-fdr_", as.character(pvalue), "-FC_",as.character(FC), "-GeneName.txt", sep="" )
      filter.name.full <- paste(as.character(edgeR.res.path), "/edgeR-downreg-", "cpmCutoff_", as.character(cpmCutoff), "-sampCutoff_", as.character(sampCutoff), "-fdr_", as.character(pvalue), "-FC_",as.character(FC), "-full.txt", sep="" )
    }
  } else if (gfilter == "both") {
    upreg <- rbind(subset(tp$table, filter==1),subset(tp$table, filter==-1))
    if (saveRes==TRUE) {
      filter.name <- paste(as.character(edgeR.res.path), "/edgeR-DEall-", "cpmCutoff_", as.character(cpmCutoff), "-sampCutoff_", as.character(sampCutoff), "-fdr_", as.character(pvalue), "-FC_",as.character(FC), "-GeneName.txt", sep="" )
      filter.name.full <- paste(as.character(edgeR.res.path), "/edgeR-DEall-", "cpmCutoff_", as.character(cpmCutoff), "-sampCutoff_", as.character(sampCutoff), "-fdr_", as.character(pvalue), "-FC_",as.character(FC), "-full.txt", sep="" )
    }
  }
  if (saveRes==TRUE) {
    if (!file.exists(as.character(filter.name))) write.table(rownames(upreg), as.character(filter.name), row.names=F, col.names=F, quote=F)
    if (!file.exists(as.character(filter.name.full))) write.table(upreg, as.character(filter.name.full), row.names=T, col.names=NA, quote=F, sep="\t")
  }
  
  return( list(filter=filter, res=tpsave$table, de=upreg, degName=rownames(upreg)) )
}

####################################################
####limma-voom DE analysis
####################################################

limmavoom.de.analysis <- function(countTabpath, metaTabpath, cpmCutoff, sampCutoff, group1, group2, padjust, Pcutoff, FCcutoff, GeneRegFilter, saveRes) {
  count.all <- read.delim(as.character(countTabpath), row.names=1, check.names=FALSE)
  metaTab <- read.delim(as.character(metaTabpath), header = T)
  metaTab <- metaTab[match(colnames(count.all), metaTab$sample),]
  
  if (dim(metaTab)[2]>2) {
    groupinfo <- metaTab[,2]
    for (i in 3:length(metaTab[1,])) {
      groupinfo <- paste(groupinfo, metaTab[,i], sep=".")
    }
    Group <- factor(groupinfo)
  } else {
    Group <- factor(metaTab[,2])
  }
  dge.count <- DGEList(counts=count.all, group=Group)
  dge.count <- calcNormFactors(dge.count)
  
  ### Remove low expression reads from NA/0 values
  cpm.count <- cpm(dge.count$counts) #cpm normalization
  threshold.value <- as.numeric(sampCutoff)
  keep = rowSums(cpm.count > as.numeric(cpmCutoff)) >=threshold.value  
  
  ##(dge.count.rmlow): DGEList of count rm low expression values.
  dge.count.rmlow <- dge.count[keep,]
  dim(dge.count.rmlow$counts)
  dge.count.rmlow$samples$lib.size <- colSums(dge.count.rmlow$counts)
  
  dge.count.rmlow <- calcNormFactors(dge.count.rmlow,method="TMM")
  log.cpm.rmlow <- cpm(dge.count.rmlow, normalized.lib.sizes=T, log=T)
  
  if (saveRes==TRUE) {
    limmavoom.res.path <- paste(getwd(), "/limma-voom-res-comp-", as.character(trim(group1)), "_", as.character(trim(group2)), sep="")
    if (!dir.exists(limmavoom.res.path)) dir.create(limmavoom.res.path)
    logcpm.fname <- paste(as.character(getwd()), "/rmlow-logcpm-", "cpmCutoof_", as.character(cpmCutoff), "-sampCutoff_", as.character(sampCutoff), "-", as.character(dim(log.cpm.rmlow)[1]), ".txt", sep=""  )
    if (!file.exists(logcpm.fname)) write.table(log.cpm.rmlow, as.character(logcpm.fname), row.names=T, col.names=NA, quote=F, sep="\t")
  
  }
  
  pvalue=as.numeric(Pcutoff)
  FC=as.numeric(FCcutoff)
  gfilter <- as.character(GeneRegFilter)
  
  f <- factor(Group, levels=levels(Group))
  design <- model.matrix(~0+f)
  rownames(design) <- rownames(dge.count.rmlow$samples)
  colnames(design) <- levels(Group)
  
  v <- voom(dge.count.rmlow, design=design, plot=T, normalize="quantile")
  fit <- lmFit(v, design)
  
  comp <- makeContrasts(contrasts=paste(as.character(trim(group2)), as.character(trim(group1)), sep="-"), levels=design )
  
  contrast.fit <- contrasts.fit(fit, contrasts=comp)
  contrast.fit <- eBayes(contrast.fit)  
  
  tpvoom.save <- topTable(contrast.fit, number=Inf, adjust="BH", sort.by="P")
  # print(head(tpvoom.save))
  
  if (saveRes==TRUE) {
    tpvoom.fname <- paste(as.character(limmavoom.res.path), "/limma-voom-full-GLMDisp-est-", "cpmCutoof_", as.character(cpmCutoff), "-sampCutoff_", as.character(sampCutoff), ".txt", sep="" )
    if(!file.exists(tpvoom.fname)) write.table(tpvoom.save, as.character(tpvoom.fname), row.names=T, col.names=NA, quote=F, sep="\t")
  }
  
  if(as.character(padjust)==T){
    tpvoom <- topTable(contrast.fit, number=Inf, adjust="BH", sort.by="none") 
    filtervoom <- decideTests(contrast.fit, adjust.method="BH", p.value=pvalue, lfc=log2(FC))
    tpvoom$filter <- as.numeric(filtervoom)
  } else {
    tpvoom <- topTable(contrast.fit, number=Inf, adjust="none", sort.by="none") 
    filtervoom <- decideTests(contrast.fit, adjust.method="none", p.value=pvalue, lfc=log2(FC))
    tpvoom$filter <- as.numeric(filtervoom)
  }
  # print(summary(filtervoom))
  
  if (gfilter == "up") {
    upreg.voom <- subset(tpvoom, filter==1)
    if (saveRes==TRUE) {
      filter.name <- paste(as.character(limmavoom.res.path), "/LimmaVoom-upreg-", "cpmCutoff_", as.character(cpmCutoff), "-sampCutoff_", as.character(sampCutoff), "-fdr_", as.character(pvalue), "-FC_",as.character(FC), "-GeneName.txt", sep="" )
      filter.name.full <- paste(as.character(limmavoom.res.path), "/LimmaVoom-upreg-", "cpmCutoff_", as.character(cpmCutoff), "-sampCutoff_", as.character(sampCutoff), "-fdr_", as.character(pvalue), "-FC_",as.character(FC), "-full.txt", sep="" )
    }
  } else if (gfilter == "down") {
    upreg.voom <- subset(tpvoom, filter==-1)
    if (saveRes==TRUE) {
      filter.name <- paste(as.character(limmavoom.res.path), "/LimmaVoom-downreg-", "cpmCutoff_", as.character(cpmCutoff), "-sampCutoff_", as.character(sampCutoff), "-fdr_", as.character(pvalue), "-FC_",as.character(FC), "-GeneName.txt", sep="" )
      filter.name.full <- paste(as.character(limmavoom.res.path), "/LimmaVoom-downreg-", "cpmCutoff_", as.character(cpmCutoff), "-sampCutoff_", as.character(sampCutoff), "-fdr_", as.character(pvalue), "-FC_",as.character(FC), "-full.txt", sep="" )
      
    }
  } else if (gfilter == "both") {
    upreg.voom <- rbind(subset(tpvoom, filter==1 ), subset(tpvoom, filter==-1) )
    if (saveRes==TRUE) {
      filter.name <- paste(as.character(limmavoom.res.path), "/LimmaVoom-DEall-", "cpmCutoff_", as.character(cpmCutoff), "-sampCutoff_", as.character(sampCutoff), "-fdr_", as.character(pvalue), "-FC_",as.character(FC), "-GeneName.txt", sep="" )
      filter.name.full <- paste(as.character(limmavoom.res.path), "/LimmaVoom-DEall-", "cpmCutoff_", as.character(cpmCutoff), "-sampCutoff_", as.character(sampCutoff), "-fdr_", as.character(pvalue), "-FC_",as.character(FC), "-full.txt", sep="" )
      
    }
  }
  
  if (saveRes==TRUE) {
    if (!file.exists(as.character(filter.name))) write.table(rownames(upreg.voom), as.character(filter.name), row.names=F, col.names=F, quote=F)
    if (!file.exists(as.character(filter.name.full))) write.table(upreg.voom, as.character(filter.name.full), row.names=T, col.names=NA, quote=F, sep="\t")
  }
  
  return( list(filter=filtervoom, res=tpvoom.save, de=upreg.voom, degName=rownames(upreg.voom)) )
}

############################################################
##### DESeq2 DE analysis
############################################################

deseq2.de.analysis <- function(countTabpath, metaTabpath, cpmCutoff, sampCutoff, group1, group2, padjust, Pcutoff, FCcutoff, GeneRegFilter, saveRes) {
  count.all <- read.delim(as.character(countTabpath), row.names=1, check.names=FALSE)
  metaTab <- read.delim(as.character(metaTabpath), header = T)
  metaTab <- metaTab[match(colnames(count.all), metaTab$sample),]
  
  if (dim(metaTab)[2]>2) {
    groupinfo <- metaTab[,2]
    for (i in 3:length(metaTab[1,])) {
      groupinfo <- paste(groupinfo, metaTab[,i], sep=".")
    }
    Group <- factor(groupinfo)
  } else {
    Group <- factor(metaTab[,2])
  }
  
  dge.count <- DGEList(counts=count.all, group=Group)
  dge.count <- calcNormFactors(dge.count)
  
  ### Remove low expression reads from NA/0 values
  cpm.count <- cpm(dge.count$counts) #cpm normalization
  threshold.value <- as.numeric(sampCutoff)
  keep = rowSums(cpm.count > as.numeric(cpmCutoff)) >=threshold.value  
  
  ##(dge.count.rmlow): DGEList of count rm low expression values.
  dge.count.rmlow <- dge.count[keep,]
  dim(dge.count.rmlow$counts)
  dge.count.rmlow$samples$lib.size <- colSums(dge.count.rmlow$counts)
  
  dge.count.rmlow <- calcNormFactors(dge.count.rmlow,method="TMM")
  log.cpm.rmlow <- cpm(dge.count.rmlow, normalized.lib.sizes=T, log=T)
  
  if (saveRes==TRUE) {
    deseq2.res.path <- paste(getwd(), "/DEseq2-res-comp-", as.character(trim(group1)), "_", as.character(trim(group2)), sep="")
    if (!dir.exists(deseq2.res.path)) dir.create(deseq2.res.path)
    logcpm.fname <- paste(as.character(getwd()), "/rmlow-logcpm-", "cpmCutoof_", as.character(cpmCutoff), "-sampCutoff_", as.character(sampCutoff), "-", as.character(dim(log.cpm.rmlow)[1]), ".txt", sep=""  )
    if (!file.exists(logcpm.fname)) write.table(log.cpm.rmlow, as.character(logcpm.fname), row.names=T, col.names=NA, quote=F, sep="\t")
    
  }
  
  pvalue=as.numeric(Pcutoff)
  FC=as.numeric(FCcutoff)
  gfilter <- as.character(GeneRegFilter)
  
  dds <- DESeqDataSetFromMatrix(countData=as.matrix(dge.count.rmlow$count), colData=DataFrame(Group), design=formula(~Group) )
  dds <- DESeq(dds, test="Wald") 
  
  res <- results(dds, contrast=c("Group", as.character(trim(group2)), as.character(trim(group1))), format="DataFrame")
  res <- as.data.frame(res)
  
  tpdeseq2.save <- res[order(res$padj) , ]
  # print(head(tpdeseq2.save))
  
  if (saveRes==TRUE) {
    tpdeseq2.fname <- paste(as.character(deseq2.res.path), "/DEseq2-full-GLMDisp-est-", "cpmCutoof_", as.character(cpmCutoff), "-sampCutoff_", as.character(sampCutoff), ".txt", sep="" )
    if(!file.exists(tpdeseq2.fname)) write.table(tpdeseq2.save, as.character(tpdeseq2.fname), row.names=T, col.names=NA, quote=F, sep="\t")
  }
  
  tpdeseq2.save$filter <- 0
  if(as.character(padjust)==T){
    tpdeseq2.save$filter[tpdeseq2.save$padj < pvalue & tpdeseq2.save$log2FoldChange > log2(FC)] <- 1
    tpdeseq2.save$filter[tpdeseq2.save$padj < pvalue & tpdeseq2.save$log2FoldChange < -log2(FC)] <- -1
  } else {
    tpdeseq2.save$filter[tpdeseq2.save$pvalue < pvalue & tpdeseq2.save$log2FoldChange > log2(FC)] <- 1
    tpdeseq2.save$filter[tpdeseq2.save$pvalue < pvalue & tpdeseq2.save$log2FoldChange < -log2(FC)] <- -1
  }
  tpdeseq2.save$filter <- as.factor(tpdeseq2.save$filter)
  # print(summary(tpdeseq2.save$filter))
  
  if (gfilter == "up") {
    upreg.deseq2 <- subset(tpdeseq2.save, filter==1)
    if (saveRes==TRUE) {
      filter.name <- paste(as.character(deseq2.res.path), "/DEseq2-upreg-", "cpmCutoff_", as.character(cpmCutoff), "-sampCutoff_", as.character(sampCutoff), "-fdr_", as.character(pvalue), "-FC_",as.character(FC), "-GeneName.txt", sep="" )
      filter.name.full <- paste(as.character(deseq2.res.path), "/DEseq2-upreg-", "cpmCutoff_", as.character(cpmCutoff), "-sampCutoff_", as.character(sampCutoff), "-fdr_", as.character(pvalue), "-FC_",as.character(FC), "-full.txt", sep="" )
    }
  } else if (gfilter == "down") {
    upreg.deseq2 <- subset(tpdeseq2.save, filter==-1)
    if (saveRes==TRUE) {
      filter.name <- paste(as.character(deseq2.res.path), "/DEseq2-downreg-", "cpmCutoff_", as.character(cpmCutoff), "-sampCutoff_", as.character(sampCutoff), "-fdr_", as.character(pvalue), "-FC_",as.character(FC), "-GeneName.txt", sep="" )
      filter.name.full <- paste(as.character(deseq2.res.path), "/DEseq2-downreg-", "cpmCutoff_", as.character(cpmCutoff), "-sampCutoff_", as.character(sampCutoff), "-fdr_", as.character(pvalue), "-FC_",as.character(FC), "-full.txt", sep="" )
    }
  } else if (gfilter == "both") {
    upreg.deseq2 <- rbind(subset(tpdeseq2.save, filter==1 ), subset(tpdeseq2.save, filter==-1) )
    if (saveRes==TRUE) {
      filter.name <- paste(as.character(deseq2.res.path), "/DEseq2-DEall-", "cpmCutoff_", as.character(cpmCutoff), "-sampCutoff_", as.character(sampCutoff), "-fdr_", as.character(pvalue), "-FC_",as.character(FC), "-GeneName.txt", sep="" )
      filter.name.full <- paste(as.character(deseq2.res.path), "/DEseq2-DEall-", "cpmCutoff_", as.character(cpmCutoff), "-sampCutoff_", as.character(sampCutoff), "-fdr_", as.character(pvalue), "-FC_",as.character(FC), "-full.txt", sep="" )
    }
  }
  if (saveRes==TRUE) {
    if (!file.exists(as.character(filter.name))) write.table(rownames(upreg.deseq2), as.character(filter.name), row.names=F, col.names=F, quote=F)
    if (!file.exists(as.character(filter.name.full))) write.table(upreg.deseq2, as.character(filter.name.full), row.names=T, col.names=NA, quote=F, sep="\t")
  }
  
  return( list(filter=(tpdeseq2.save$filter), res=tpdeseq2.save, de=upreg.deseq2, degName=rownames(upreg.deseq2)) )
}

###############################################################
##Cuffdiff res process
###############################################################

cuffdiff.res.process <- function(respath, group1, group2, padjust, Pcutoff, FCcutoff, GeneRegFilter, saveRes) {
  cuffdiff.res<- read.table(as.character(respath), header = T)
  head(cuffdiff.res)
  
  if(saveRes==TRUE){
    cuff.res.path <- paste(getwd(), "/cuffdiff-res", sep="")
    if (!dir.exists(cuff.res.path)) dir.create(cuff.res.path)
  }
  pvalue=as.numeric(Pcutoff)
  FC=as.numeric(FCcutoff)
  gfilter <- as.character(GeneRegFilter)
  
  cuffdiff.res.sub <- cuffdiff.res[ ( (cuffdiff.res$sample_1 == as.character(group1) & cuffdiff.res$sample_2==as.character(group2)) | (cuffdiff.res$sample_1 == as.character(group2) & cuffdiff.res$sample_2==as.character(group1))), ]
  cuffdiff.res.sub.order <- cuffdiff.res.sub[ order(cuffdiff.res.sub$q_value),]
  
  if(saveRes==TRUE) {
    cuff.save.name <- paste(as.character(cuff.res.path), "/cuffdiff-", as.character(group1), "_", as.character(group2), ".txt", sep="")
    if (!file.exists(cuff.save.name)) write.table(cuffdiff.res.sub.order, as.character(cuff.save.name), row.names = F, quote=F, sep="\t")
  }
  
  cuffdiff.res.sub.order$filter <- 0
  if(as.character(padjust)==T){
    cuffdiff.res.sub.order$filter[cuffdiff.res.sub.order$q_value < pvalue & cuffdiff.res.sub.order$log2.fold_change. > log2(FC)] <- 1
    cuffdiff.res.sub.order$filter[cuffdiff.res.sub.order$q_value < pvalue & cuffdiff.res.sub.order$log2.fold_change. < -log2(FC)] <- -1
  } else {
    cuffdiff.res.sub.order$filter[cuffdiff.res.sub.order$q_value < pvalue & cuffdiff.res.sub.order$log2.fold_change. > log2(FC)] <- 1
    cuffdiff.res.sub.order$filter[cuffdiff.res.sub.order$q_value < pvalue & cuffdiff.res.sub.order$log2.fold_change. < -log2(FC)] <- -1
  }
  cuffdiff.res.sub.order$filter <- as.factor(cuffdiff.res.sub.order$filter)
  # print(summary(cuffdiff.res.sub.order$filter))
  
  if (gfilter == "up") {
    upreg.cuffdiff <- subset(cuffdiff.res.sub.order, filter==1)
    if (saveRes==TRUE) {
      filter.name <- paste(as.character(cuff.res.path), "/cuffdiff-upreg-", "fdr_", as.character(pvalue), "-FC_",as.character(FC), "-GeneName.txt", sep="" )
      filter.name.full <- paste(as.character(cuff.res.path), "/cuffdiff-upreg-", "fdr_", as.character(pvalue), "-FC_",as.character(FC), "-full.txt", sep="" )
    }
  } else if (gfilter == "down") {
    upreg.cuffdiff <- subset(cuffdiff.res.sub.order, filter==-1)
    if (saveRes==TRUE) {
      filter.name <- paste(as.character(cuff.res.path), "/cuffdiff-downreg-", "fdr_", as.character(pvalue), "-FC_",as.character(FC), "-GeneName.txt", sep="" )
      filter.name.full <- paste(as.character(cuff.res.path), "/cuffdiff-downreg-", "fdr_", as.character(pvalue), "-FC_",as.character(FC), "-full.txt", sep="" )
    }
  } else if (gfilter == "both") {
    upreg.cuffdiff <- rbind(subset(cuffdiff.res.sub.order, filter==1 ), subset(cuffdiff.res.sub.order, filter==-1) )
    if (saveRes==TRUE) {
      filter.name <- paste(as.character(cuff.res.path), "/cuffdiff-DEall-", "fdr_", as.character(pvalue), "-FC_",as.character(FC), "-GeneName.txt", sep="" )
      filter.name.full <- paste(as.character(cuff.res.path), "/cuffdiff-DEall-", "fdr_", as.character(pvalue), "-FC_",as.character(FC), "-full.txt", sep="" )
    }
  }
  if (saveRes==TRUE) {
    if (!file.exists(as.character(filter.name))) write.table(upreg.cuffdiff$gene, as.character(filter.name), row.names=F, col.names=F, quote=F)
    if (!file.exists(as.character(filter.name.full))) write.table(upreg.cuffdiff, as.character(filter.name.full), row.names=F, quote=F, sep="\t")
  }
  
  return( list(filter=(cuffdiff.res.sub.order$filter), res=cuffdiff.res.sub.order, de=upreg.cuffdiff, degName=as.character(upreg.cuffdiff$gene)) )
  
}
