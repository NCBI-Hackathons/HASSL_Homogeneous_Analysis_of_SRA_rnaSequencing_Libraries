#Load packages
library(DESeq2)
library("ggplot2")

#Read in data
countData <- read.table("/Users/dillmana/Documents/courses/merge_genecount.txt", header=TRUE, stringsAsFactors=TRUE, row.names="Name")

#filter out low coverage reads
rs <- rowSums ( countData )
use <- (rs > quantile(rs, 0.4))
table(use)
cdsFilt <- countData[ use, ]

#name of files
colData <- read.csv("Design.csv", header=TRUE, stringsAsFactors=TRUE, row.names="sample_ID")

#combine tables and normalize for library depth
ddsfil <- DESeqDataSetFromMatrix(countData = cdsFilt, colData = colData, design = ~ Sample)
ddsfil <- estimateSizeFactors(ddsfil)
deseq_ncounts <- counts(ddsfil, normalized=TRUE)
genes <- row.names(deseq_ncounts)
tdata <- as.data.frame(t(deseq_ncounts))
colnames(tdata) <- genes

#make plots for all genes
plotHistFunc <- function(x, na.rm = TRUE, ...) {
  nm <- names(x)
  for (i in seq_along(nm)) {
    plots <-ggplot(x,aes_string(x = nm[i])) + geom_density(alpha = .5,fill = "dodgerblue")
    ggsave(plots,filename=paste(nm[i],".jpeg",sep=""))
  }
}
plotHistFunc(tdata)


