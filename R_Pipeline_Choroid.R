library("BiocManager")
library("DESeq2")
library(reshape2)
library(ggplot2)
library("RColorBrewer") # Load a package giving more colors
library("pheatmap") # load a package for making heatmaps

#Read the data into Data frame
setwd("/storage/hpc/data/bourasm/CXCR5")
CXCR5<-read.table("C57_vs_CXCR5_Choroid_results.txt", header = TRUE, row.names = 1)
CXCR5Design <-data.frame(
  row.names = colnames(CXCR5),
  condition = c("untraited", "untraited", "untraited",
                "traited", "traited", "traited"),
  libType = c("paired-end","paired-end","paired-end",
              "paired-end","paired-end","paired-end"))

pf=function (x,y){
  pdf(x) 
  y 
  dev.off()
}

pf("results_Boxplot_CXCR5.pdf", boxplot(CXCR5) )
pf("results_Hist_CXCR5.pdf", hist(CXCR5[,1]) )# Plotting only the first sample (column 1)


pseudoCount = log2(CXCR5 + 1) # log-transform to make numbers on scale (+1 to avoid zeroes)
boxplot(pseudoCount)
pf("results_Boxplot_LogCount.pdf", boxplot(pseudoCount))

pf("results_Hist_LogCount.pdf", hist(pseudoCount[,1]))

pseudoCount = as.data.frame(pseudoCount)
df = melt(pseudoCount, variable.name = "Samples", value.name = "count") # reshape the matrix 
df01 = data.frame(df, Condition = substr(df$Samples, 1, 6))

p<- ggplot(df01, aes(x = df$Samples, y = count, fill = Condition)) + geom_boxplot() + xlab("") +
  ylab(expression(log[2](count + 1)))

pf("results_Boxplot_LogCount_02.pdf", print(p))



p<-ggplot(df01, aes(x = count, colour = Samples, fill = Samples)) + ylim(c(0, 0.17)) +
  geom_density(alpha = 0.2, size = 1.25) + facet_wrap(~ Condition) +
  theme(legend.position = "top") + xlab(expression(log[2](count + 1)))

pf("results_Boxplot_LogCount_03.pdf", print(p))

dds = DESeqDataSetFromMatrix(countData = CXCR5,
                             colData = CXCR5Design,
                             design = ~ condition) # we're testing for the different condidtions
dds

rld = rlogTransformation(dds)

#Need to create a more general function to plotting.
pdf("results_Boxplot_DESeqLog.pdf") 
   plot(log2( 1 + counts(dds)[ , 1:2] ),
     pch=16, cex=0.3, main = "log2")
dev.off()

pdf("results_Boxplot_DESeqrld01.pdf")
plot(assay(rld)[ , 1:2], # The assay function returns the count values of rld in a matrix
     pch=16, cex=0.3, main = "rlog")
dev.off()

distsRL <- dist(t(assay(rld))) # Calculate distances using transformed (and normalized) counts
mat <- as.matrix(distsRL) # convert to matrix
rownames(mat) <- colnames(mat) <- with(colData(dds), paste(condition)) # set rownames in the matrix
colnames(mat) = NULL # remove column names
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) # set colors
pdf("pheatmap.pdf")
   pheatmap(mat,
         clustering_distance_rows=distsRL,
         clustering_distance_cols=distsRL,
         col=colors)
dev.off()

pdf("plotPCArld.pdf")
plotPCA(rld, intgroup=c("condition"))
dev.off()

pdf("plotpseudoCount_CX5_C2.pdf")
plot(pseudoCount[,2], pseudoCount[,5], pch = 20, xlab = "C57_C2", ylab = "CX5_C2") + # pch specifies the type of symbols. Filled dots in this case.
  abline(0,1, col = "red") # make a diagonal line
dev.off()

x = pseudoCount[, 2] # extract fourth column
y = pseudoCount[, 5] # extract seventh column
M = x - y # M-values (differences between samples)
A = (x + y)/2 # A-values (averages)
df = data.frame(A, M)

p<-ggplot(df, aes(x = A, y = M)) + geom_point(size = 1.5, alpha = 1/5) +
  geom_hline(aes(yintercept = 0), color = "blue3") + stat_smooth(se =
  FALSE, method = "loess", color = "red3") + 
  xlab("Average gene count") +
  ylab("Count difference") +
  ggtitle("C57_C2 vs. CX5_C2")

pf("AverageGeneCount.pdf",p)


dds = estimateSizeFactors(dds)
sizeFactors(dds)


norm_counts = counts(dds, normalized = TRUE) # Extract the normalized counts
pseudoCount = log2(norm_counts + 1) # convert to log-scale for visualization
df = melt(pseudoCount) # transpose the matrix
df = data.frame(df, Condition = substr(df$Var2, 1, 4))

pdf("AverageGeneCountLog.pdf")
ggplot(df, aes(x = df$Var2, y = value, fill = Condition)) + geom_boxplot() + xlab("") +
  ylab(expression(log[2](count + 1)))
dev.off()

dds = DESeq(dds)
res <- results(dds)
mcols(res, use.names=TRUE)
summary(res)

pdf("DESeq_dds.pdf")
plotMA(res, ylim=c(-7,7))
dev.off()

resShrink = lfcShrink(dds, coef=2)

pdf("DESeq_dds_ver02.pdf")
plotMA(resShrink, ylim=c(-5,5))
dev.off()

resSig <- subset(res, padj < 0.1)
resSig

head(resShrink[order(resShrink$log2FoldChange), ], 10)

attributes(resSig)[4]


View(head(resShrink[ order(resShrink$log2FoldChange), ], 10))
pdf("ENSMUSG00000107651.pdf")
plotCounts(dds, "ENSMUSG00000107651", "condition")

head(resShrink[ order(resShrink$log2FoldChange, decreasing=TRUE), ], 10)
pdf("ENSMUSG00000020866.pdf")
plotCounts(dds, "ENSMUSG00000020866", "condition")

head(res[order(res$padj),], 5) # order by padjusted and print the top 5
pdf("ENSMUSG00000059412.pdf")
plotCounts(dds, "ENSMUSG00000059412", "condition")



resLFC1 = results(dds, lfcThreshold = 1)
summary(resLFC1)

pdf("resLFC1.pdf")
plotMA(resLFC1, ylim=c(-7,7)) +
  abline(h = 1, col = "blue") +
  abline(h = -1, col = "blue")


resLFC1Srhunk = lfcShrink(dds, coef=2, res=resLFC1)

pdf("resLFC1Srhunk.pdf")
plotMA(resLFC1Srhunk, ylim=c(-5,5))+
  abline(h = 1, col = "blue") +
  abline(h = -1, col = "blue")

res.05 <- results(dds, alpha=.05)
table(res.05$padj < .05)
summary(res.05)

res.05Shrunk = lfcShrink(dds, coef = 2, res=res.05)
pdf("res.05Shrunk.pdf")
plotMA(res.05Shrunk, alpha = 0.05, ylim=c(-5,5)) +
  abline(h = 1, col = "blue") +
  abline(h = -1, col = "blue")

resUpReg = res[which(res$log2FoldChange < 0), ] # get the upregulated genes
head(resUpReg[order(resUpReg$padj),], 5) # order by padjusted and print the top 5

pdf("ENSMUSG00000072476.pdf")
plotCounts(dds, "ENSMUSG00000072476", "condition")

mat = assay(rld)[ head(order(res$padj),30), ] # select the top 30 genes with the lowest padj
mat = mat - rowMeans(mat) # Subtract the row means from each value
# Optional, but to make the plot nicer:
df = as.data.frame(colData(rld)[,c("condition")]) # Create a dataframe with a column of the conditions
colnames(df) = "condition" # Rename the column header
rownames(df) = colnames(mat) # add rownames
# and plot the actual heatmap
pdf("pheatmapDESeq_ver02.pdf")
pheatmap(mat, annotation_col=df)
dev.off()

assay(rld)["ENSMUSG00000072476",]
