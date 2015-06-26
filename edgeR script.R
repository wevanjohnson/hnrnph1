expr <- read.table("htseq_exprMatrix.txt", header=T)

congenic <- rep(c("no", "yes"), each=2, length=32)
sample <- rep(1:16, each=2)
cage <- rep(c(89,87,81,80,89,87,86,90,90,85,90,84,83,82,83,82), each=2)
lane <- c(rep(c(5,6),8), rep(c(7,8),8))

pheno <- data.frame(congenic=congenic, sample=sample, cage=cage, lane=lane)
pheno <- pheno[-25,]
row.names(pheno) <- names(expr)

library(edgeR)
keep <- rowSums(cpm(expr)>1) >= 16
expr <- expr[keep, ] #dim [1] 14177    31
table(keep)

y <- DGEList(counts=expr, group=as.factor(pheno$congenic))
y <- calcNormFactors(y)

# No adjustment for cage effect
design <- model.matrix(~as.factor(pheno$congenic))
rownames(design) <- colnames(y)
y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
#plotBCV(y1)
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
topGenes <- topTags(lrt,n=nrow(y))

# fold-change > 1.5 or < 1.5 & FDR < 0.05
topGenes_FC1.5_fdr <- subset(topGenes$table, (logFC < log(1/1.5, base=2) | logFC > log(1.5, base=2)) & FDR < 0.05)

# FDR < 0.05
topGenes_fdr <- subset(topGenes$table, FDR < 0.05)

# FDR < 0.01
topGenes_fdr0.01 <- subset(topGenes$table, FDR < 0.01)

#plots
summary(dt <- decideTestsDGE(lrt))
isDE <- as.logical(dt)
DEnames <- rownames(y)[isDE]
plotSmear(lrt, de.tags=DEnames)
abline(h=c(-1,1), col="blue")

# adjustment for cage effect
design <- model.matrix(~as.factor(pheno$cage)+as.factor(pheno$congenic))
rownames(design) <- colnames(y)
y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
#plotBCV(y)
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
topGenes_cage <- topTags(lrt,n=nrow(y))

# fold-change > 1.5 or < 1.5 & FDR < 0.05
topGenes_cage_FC1.5_fdr <- subset(topGenes_cage$table, (logFC < log(1/1.5, base=2) | logFC > log(1.5, base=2)) & FDR < 0.05)

# FDR < 0.05
topGenes_cage_fdr <- subset(topGenes_cage$table, FDR < 0.05)

# FDR < 0.01
topGenes_cage_fdr0.01 <- subset(topGenes_cage$table, FDR < 0.01)

# FDR < 0.1
topGenes_cage_fdr0.1 <- subset(topGenes_cage$table, FDR < 0.1)

# FDR < 0.2
topGenes_cage_fdr0.2 <- subset(topGenes_cage$table, FDR < 0.2)


#plots
summary(dt <- decideTestsDGE(lrt))
isDE <- as.logical(dt)
DEnames <- rownames(y)[isDE]
plotSmear(lrt, de.tags=DEnames)
abline(h=c(-1,1), col="blue")

# output results
write.csv(topGenes_FC1.5_fdr, file="noAdjCage_FDR0.05_FC1.5.csv")
write.csv(topGenes_fdr, file="noAdjCage_FDR0.05.csv")
write.csv(topGenes_fdr0.01, file="noAdjCage_FDR0.01.csv")

write.csv(topGenes_cage_FC1.5_fdr, file="AdjCage_FDR0.05_FC1.5.csv")
write.csv(topGenes_cage_fdr, file="AdjCage_FDR0.05.csv")
write.csv(topGenes_cage_fdr0.01, file="AdjCage_FDR0.01.csv")
write.csv(topGenes_cage_fdr0.1, file="AdjCage_FDR0.1.csv")
write.csv(topGenes_cage_fdr0.2, file="AdjCage_FDR0.2.csv")





