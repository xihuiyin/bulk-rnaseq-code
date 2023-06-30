library("DESeq2")

directory <- "../htseq/htseq_results"
sampleFiles <- c('lisa_1UTX.txt','lisa_2WT.txt','lisa_5UTX.txt','lisa_6WT.txt','lisa_8WT.txt','lisa_9UTX.txt','lisa_10WT.txt')
condition <- c("UTX","WT","UTX","WT","WT","UTX","WT")
sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = condition)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)

keep <- rowSums(counts(ddsHTSeq)) >= 10
ddsHTSeq <- ddsHTSeq[keep,]


ddsHTSeq <- DESeq(ddsHTSeq)
res <- results(ddsHTSeq)
res


library(apeglm)
res <- results(ddsHTSeq, contrast=c("condition","WT","UTX"))

resLFC <- lfcShrink(ddsHTSeq, coef="condition_WT_vs_UTX", type="apeglm")
resLFC

resOrdered <- res[order(res$pvalue),]
summary(res)

plotMA(resLFC, ylim=c(-2,2))
idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]

# because we are interested in treated vs untreated, we set 'coef=2'
resNorm <- lfcShrink(ddsHTSeq, coef=2, type="normal")
resAsh <- lfcShrink(ddsHTSeq, coef=2, type="ashr")

par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")

vsd <- vst(ddsHTSeq, blind=FALSE)
plotPCA(vsd, intgroup=c("condition"))
resLFC <- resLFC[order(resLFC$padj),]
head(resLFC)

par(mfrow=c(2,3))

plotCounts(ddsHTSeq, gene="ENSMUSG00000074800", intgroup="condition")
plotCounts(ddsHTSeq, gene="ENSMUSG00000032691", intgroup="condition")
plotCounts(ddsHTSeq, gene="ENSMUSG00000044550", intgroup="condition")
plotCounts(ddsHTSeq, gene="ENSMUSG00000072717", intgroup="condition")
plotCounts(ddsHTSeq, gene="ENSMUSG00000104915", intgroup="condition")
plotCounts(ddsHTSeq, gene="ENSMUSG00000079853", intgroup="condition")

plotCounts(ddsHTSeq, gene="ENSMUSG00000094087", intgroup="condition")
plotCounts(ddsHTSeq, gene="ENSMUSG00000002012", intgroup="condition")
plotCounts(ddsHTSeq, gene="ENSMUSG00000080727", intgroup="condition")
plotCounts(ddsHTSeq, gene="ENSMUSG00000039114", intgroup="condition")
plotCounts(ddsHTSeq, gene="ENSMUSG00000025480", intgroup="condition")

#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(resLFC, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(resLFC, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(resLFC, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
