####volcano plot function
volcano.plot.fn.deseq2 <- function(input.fname, fname, pval, fc) {
  plotres <- read.delim(input.fname, header = T, row.names = 1, check.names = F)
  plotres$regulation <- 0
  plotres$regulation[plotres$padj<=pval & plotres$log2FoldChange >= log2(fc) ] <- 1
  plotres$regulation[plotres$padj<=pval & plotres$log2FoldChange <= -log2(fc) ] <- -1
  print("DESeq2 Res: ")
  print(table(plotres$regulation))
  print("======")
  # plotres$padj[is.na(plotres$padj)] <- 1 ##log2 FC of these DEGs are large
  xlim.max <- max(-min(plotres$log2FoldChange), max(plotres$log2FoldChange))
  ylim.max <- max(-log10(plotres$padj+1e-10), na.rm = T) + 0.5
  library(ggplot2)
  g <- ggplot(data=plotres, aes(x=log2FoldChange, y=-log10(padj+1e-10), colour=factor(regulation))) + geom_point()
  g <- g + labs(x="log2 FC", y="-log10 (FDR)")
  g <- g + scale_size_manual(values=c("0"=0.5, "1"=2, "-1"=2 ))
  g <- g + xlim(-xlim.max, xlim.max) 
  g <- g + geom_hline(yintercept=-log10(pval)) 
  g <- g + geom_vline(xintercept=log2(fc)) + geom_vline(xintercept=-log2(fc))
  g <- g + theme(legend.title=element_blank())
  g <- g + scale_color_manual(values = c("0"="grey70", "1"="#EE7621", "-1"="#191970" ), 
                              breaks=c("1", "-1", "0"),
                              labels=c("Up-regulated", "Down-regulated", "non DE"))
  # g <- g + theme(legend.title=element_blank(),legend.position = "top", legend.direction="vertical")
  g <- g + theme(axis.title.x=element_text(size=20)) + theme(axis.title.y=element_text(size=20,vjust=1.5))
  g <- g + theme(axis.text.x=element_text(size=20)) + theme(axis.text.y=element_text(size=20))
  g <- g + theme(legend.text=element_text(size=20))
  g <- g + guides(colour = guide_legend(override.aes = list(size=10)))
  g <- g + annotate("text",  x=-xlim.max+0.2, y = ylim.max,  
                    label = as.character(sum(plotres$regulation==-1)),
                    size = 9, color = "#191970")
  g <- g + annotate("text",  x= xlim.max-0.2, y = ylim.max, 
                    label = as.character(sum(plotres$regulation==1)),
                    size = 9, color = "#EE7621")
  ggsave(filename = as.character(fname), plot = g, height = 5, width = 7)
}

volcano.plot.fn.edger <- function(input.fname, fname, pval, fc) {
  plotres <- read.delim(input.fname, header = T, row.names = 1, check.names = F)
  plotres$regulation <- 0
  plotres$regulation[plotres$FDR<=pval & plotres$logFC >= log2(fc) ] <- 1
  plotres$regulation[plotres$FDR<=pval & plotres$logFC <= -log2(fc) ] <- -1
  print("edgeR Res: ")
  print(table(plotres$regulation))
  print("======")
  # plotres$FDR[is.na(plotres$FDR)] <- 1 ##log2 FC of these DEGs are large
  xlim.max <- max(-min(plotres$logFC), max(plotres$logFC))
  ylim.max <- max(-log10(plotres$FDR+1e-10), na.rm = T) + 0.5
  library(ggplot2)
  g <- ggplot(data=plotres, aes(x=logFC, y=-log10(FDR+1e-10), colour=factor(regulation))) + geom_point()
  g <- g + labs(x="log2 FC", y="-log10 (FDR)")
  g <- g + scale_size_manual(values=c("0"=0.5, "1"=2, "-1"=2 ))
  g <- g + xlim(-xlim.max, xlim.max) 
  g <- g + geom_hline(yintercept=-log10(pval)) 
  g <- g + geom_vline(xintercept=log2(fc)) + geom_vline(xintercept=-log2(fc))
  g <- g + theme(legend.title=element_blank())
  g <- g + scale_color_manual(values = c("0"="grey70", "1"="#EE7621", "-1"="#191970" ), 
                              breaks=c("1", "-1", "0"),
                              labels=c("Up-regulated", "Down-regulated", "non DE"))
  # g <- g + theme(legend.title=element_blank(),legend.position = "top", legend.direction="vertical")
  g <- g + theme(axis.title.x=element_text(size=20)) + theme(axis.title.y=element_text(size=20,vjust=1.5))
  g <- g + theme(axis.text.x=element_text(size=20)) + theme(axis.text.y=element_text(size=20))
  g <- g + theme(legend.text=element_text(size=20))
  g <- g + guides(colour = guide_legend(override.aes = list(size=10)))
  g <- g + annotate("text",  x=-xlim.max+1, y = ylim.max,  
                    label = as.character(sum(plotres$regulation==-1)),
                    size = 9, color = "#191970")
  g <- g + annotate("text",  x= xlim.max-1, y = ylim.max, 
                    label = as.character(sum(plotres$regulation==1)),
                    size = 9, color = "#EE7621")
  ggsave(filename = as.character(fname), plot = g, height = 5, width = 7)
}

volcano.plot.fn.voom <- function(input.fname, fname, pval, fc) {
  plotres <- read.delim(input.fname, header = T, row.names = 1, check.names = F)
  plotres$regulation <- 0
  plotres$regulation[plotres$adj.P.Val<=pval & plotres$logFC >= log2(fc) ] <- 1
  plotres$regulation[plotres$adj.P.Val<=pval & plotres$logFC <= -log2(fc) ] <- -1
  print("Limma voom Res: ")
  print(table(plotres$regulation))
  print("======")
  # plotres$adj.P.Val[is.na(plotres$adj.P.Val)] <- 1 ##log2 FC of these DEGs are large
  xlim.max <- max(-min(plotres$logFC), max(plotres$logFC))
  ylim.max <- max(-log10(plotres$adj.P.Val+1e-10), na.rm = T) + 0.5
  library(ggplot2)
  g <- ggplot(data=plotres, aes(x=logFC, y=-log10(adj.P.Val+1e-10), colour=factor(regulation))) + geom_point()
  g <- g + labs(x="log2 FC", y="-log10 (FDR)")
  g <- g + scale_size_manual(values=c("0"=0.5, "1"=2, "-1"=2 ))
  g <- g + xlim(-xlim.max, xlim.max) 
  g <- g + geom_hline(yintercept=-log10(pval)) 
  g <- g + geom_vline(xintercept=log2(fc)) + geom_vline(xintercept=-log2(fc))
  g <- g + theme(legend.title=element_blank())
  g <- g + scale_color_manual(values = c("0"="grey70", "1"="#EE7621", "-1"="#191970" ), 
                              breaks=c("1", "-1", "0"),
                              labels=c("Up-regulated", "Down-regulated", "non DE"))
  # g <- g + theme(legend.title=element_blank(),legend.position = "top", legend.direction="vertical")
  g <- g + theme(axis.title.x=element_text(size=20)) + theme(axis.title.y=element_text(size=20,vjust=1.5))
  g <- g + theme(axis.text.x=element_text(size=20)) + theme(axis.text.y=element_text(size=20))
  g <- g + theme(legend.text=element_text(size=20))
  g <- g + guides(colour = guide_legend(override.aes = list(size=10)))
  g <- g + annotate("text",  x=-xlim.max+1, y = ylim.max,  
                    label = as.character(sum(plotres$regulation==-1)),
                    size = 9, color = "#191970")
  g <- g + annotate("text",  x= xlim.max-1, y = ylim.max, 
                    label = as.character(sum(plotres$regulation==1)),
                    size = 9, color = "#EE7621")
  ggsave(filename = as.character(fname), plot = g, height = 5, width = 7)
}

get.deg <- function(input.fname, fc, pval, method) {
  res <- read.delim(input.fname, header = T, row.names = 1, check.names = F)
  res$filter <- 0
  if (method == 'edgeR') {
    res$filter[res$FDR <= pval & res$logFC >= log2(fc)] <- 1
    res$filter[res$FDR <= pval & res$logFC <= -log2(fc)] <- -1
  } else if (method == 'DESeq2') {
    res$filter[res$padj <= pval & res$log2FoldChange >= log2(fc)] <- 1
    res$filter[res$padj <= pval & res$log2FoldChange <= -log2(fc)] <- -1
  } else if (method == 'voom') {
    res$filter[res$adj.P.Val <= pval & res$logFC >= log2(fc)] <- 1
    res$filter[res$adj.P.Val <= pval & res$logFC <= -log2(fc)] <- -1
  }
  res.output <- res[res$filter!=0,]
  return(res.output)
}
####End functions setting up
####################





