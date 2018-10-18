library(EBSeqHMM)
library(EBSeq)
library(edgeR)
library(sva)
library(biomaRt)
library(RColorBrewer)
library(limma)
library(NOISeq)

## load & prepare dataset ##
aitordata <- read.table("/media/shinigami/Elements/RNAseq/countsaitor.txt", comment.char="#")
aitordata <- as.matrix(aitordata)
rownames(aitordata) <- aitordata[,1]
aitordata <- aitordata[,-1]
aitordata <- aitordata[-1,]
exprMatrix <- aitordata[,6:29]
colnames(exprMatrix) <- c("D8R1", "D24R3", "D24R2", "D8R3", "D8R2", "D24R1", "M24R3", "M8R2", "M8R1", "M8R3", "M24R2", "M24R1", "T8R1", "T24R2", "T24R3", "T24R1", "T8R2", "T8R3", "TM8R3", "TM8R1", "TM24R2", "TM24R1", "TM8R2", "TM24R3")
mode(exprMatrix) <- 'numeric'

## filtering ##
keep <- filterByExpr(exprMatrix)
exprMatrixflt <- exprMatrix[keep,]
checker <- aveLogCPM(exprMatrixflt)
expressed <- exprMatrixflt[checker>1,]

## phenotype matrix ##
samples_series <- colnames(exprMatrixflt)
replicate_series <- c('R1','R3','R2','R3','R2','R1','R3','R2','R1','R3','R2','R1','R1','R2','R3','R1','R2','R3','R3','R1','R2','R1','R2','R3')
timepoint_series <- c('8','24','24','8','8','24','24','8','8','8','24','24','8','24','24','24','8','8','8','8','24','24','8','24')
condition_series <- c("D8", "D24", "D24", "D8", "D8", "D24", "M24", "M8", "M8", "M8", "M24", "M24", "T8", "T24", "T24", "T24", "T8", "T8", "TM8", "TM8", "TM24", "TM24", "TM8", "TM24")
phenodata <- as.data.frame(cbind(samples_series,condition_series,replicate_series,timepoint_series))
rownames(phenodata) <- phenodata[,1]
phenodata <- phenodata[,-1]
colnames(phenodata)<-c('condition','replicate','timepoint')

nsamples <- ncol(expressed)
col <- brewer.pal(nsamples, "Paired")
boxplot(log2(expressed), col=col, main="")
title(main="Log2 Raw data",ylab="Log-counts")

## Normalization (Median, Quantile, TMM) ##
sizes <- MedianNorm(expressed)
GeneNormData <- GetNormalizedMat(expressed, sizes)
GeneNormData2 <- normalizeQuantiles(expressed)
GeneNormData3 <- tmm(expressed, long = 1000, lc = 0, k = 0)

boxplot(log2(GeneNormData), las=2, col=col, main="")
title(main="Median Normalised data",ylab="Log2-counts")
plotMDS(log2(GeneNormData), top = 1000, labels = NULL, col = as.numeric(phenodata$condition), cex = 2)
title(main="Treatment groups - Median")
plotMDS(log2(GeneNormData), top = 1000, labels = NULL, col = as.numeric(phenodata$replicate), cex = 2)
title(main="Replicate groups - Median")

boxplot(log2(GeneNormData2), las=2, col=col, main="")
title(main="Quantile Normalised data",ylab="Log2-counts")
plotMDS(log2(GeneNormData2), top = 1000, labels = NULL, col = as.numeric(phenodata$condition), cex = 2)
title(main="Treatment groups - Quantile")
plotMDS(log2(GeneNormData2), top = 1000, labels = NULL, col = as.numeric(phenodata$replicate), cex = 2)
title(main="Replicate groups - Quantile")

boxplot(log2(GeneNormData3), las=2, col=col, main="")
title(main="TMM Normalised data",ylab="Log2-counts")
plotMDS(log2(GeneNormData3), top = 1000, labels = NULL, col = as.numeric(phenodata$condition), cex = 2)
title(main="Treatment groups - TMM")
plotMDS(log2(GeneNormData3), top = 1000, labels = NULL, col = as.numeric(phenodata$replicate), cex = 2)
title(main="Replicate groups - TMM")

## Batch effect adjustment (Combat) ##
logCPM <- cpm(GeneNormData, log = T, prior.count = .5)
logCPM2 <- cpm(GeneNormData2, log = T, prior.count = .5)
logCPM3 <- cpm(GeneNormData3, log = T, prior.count = .5)

batch = phenodata$replicate
modcombat = model.matrix(~condition, data=phenodata)
cleaned = ComBat(dat=logCPM, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
cleaned2 = ComBat(dat=logCPM2, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
cleaned3 = ComBat(dat=logCPM3, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

boxplot(cleaned, las=2, col=col, main="")
title(main="Median+Combat data",ylab="Log2-cpm")
plotMDS(cleaned, top = 1000, labels = NULL, col = as.numeric(phenodata$condition), cex = 2)
title(main="Treatment groups - Median")
plotMDS(cleaned, top = 1000, labels = NULL, col = as.numeric(phenodata$replicate), cex = 2)
title(main="Replicate groups - Median")

boxplot(cleaned2, las=2, col=col, main="")
title(main="Quantile+Combat data",ylab="Log2-cpm")
plotMDS(cleaned2, top = 1000, labels = NULL, col = as.numeric(phenodata$condition), cex = 2)
title(main="Treatment groups - Quantile")
plotMDS(cleaned2, top = 1000, labels = NULL, col = as.numeric(phenodata$replicate), cex = 2)
title(main="Replicate groups - Quantile")

boxplot(cleaned3, las=2, col=col, main="")
title(main="TMM+Combat data",ylab="Log2-cpm")
plotMDS(cleaned3, top = 1000, labels = NULL, col = as.numeric(phenodata$condition), cex = 2)
title(main="Treatment groups - TMM")
plotMDS(cleaned3, top = 1000, labels = NULL, col = as.numeric(phenodata$replicate), cex = 2)
title(main="Replicate groups - TMM")

## data transformation ##
antilog<-function(lx,base) {
  lbx<-lx/log(exp(1),base=base)
  result<-exp(lbx)
  result
}
cleanedanti <- antilog(cleaned, 2)
exprset_8h <- cleanedanti[,phenodata$timepoint == 8]
exprset_24h <- cleanedanti[,phenodata$timepoint == 24]
cleanedanti2 <- antilog(cleaned2, 2)
exprset2_8h <- cleanedanti2[,phenodata$timepoint == 8]
exprset2_24h <- cleanedanti2[,phenodata$timepoint == 24]
cleanedanti3 <- antilog(cleaned3, 2)
exprset3_8h <- cleanedanti3[,phenodata$timepoint == 8]
exprset3_24h <- cleanedanti3[,phenodata$timepoint == 24]

## EBSEQHMM ##
sf <- rep(1,12)
conditions <- as.factor(rep(c('dmso','mkc','taxol','taxolmkc'), each=3))

EBHMMNBGeneout_8h <- EBHMMNBMultiEM_2chain(Data=exprset_8h, sizeFactors=sf, Conditions=conditions, UpdateRd=10, PriorFC=1.3)
GeneConfCalls_8h <- GetConfidentCalls(EBHMMNBGeneout_8h, FDR=.05, cutoff=.5, OnlyDynamic=T)
results_8h <- as.data.frame(GeneConfCalls_8h$Overall)
table(results_8h$Most_Likely_Path)
EBHMMNBGeneout_24h <- EBHMMNBMultiEM_2chain(Data=exprset_24h, sizeFactors=sf, Conditions=conditions, UpdateRd=10, PriorFC=1.3)
GeneConfCalls_24h <- GetConfidentCalls(EBHMMNBGeneout_24h, FDR=.05, cutoff=.5, OnlyDynamic=T)
results_24h <- as.data.frame(GeneConfCalls_24h$Overall)
table(results_24h$Most_Likely_Path)

EBHMMNBGeneout2_8h <- EBHMMNBMultiEM_2chain(Data=exprset2_8h, sizeFactors=sf, Conditions=conditions, UpdateRd=10, PriorFC=1.3)
GeneConfCalls2_8h <- GetConfidentCalls(EBHMMNBGeneout2_8h, FDR=.05, cutoff=.5, OnlyDynamic=T)
results2_8h <- as.data.frame(GeneConfCalls2_8h$Overall)
table(results2_8h$Most_Likely_Path)
EBHMMNBGeneout2_24h <- EBHMMNBMultiEM_2chain(Data=exprset2_24h, sizeFactors=sf, Conditions=conditions, UpdateRd=10, PriorFC=1.3)
GeneConfCalls2_24h <- GetConfidentCalls(EBHMMNBGeneout2_24h, FDR=.05, cutoff=.5, OnlyDynamic=T)
results2_24h <- as.data.frame(GeneConfCalls2_24h$Overall)
table(results2_24h$Most_Likely_Path)

EBHMMNBGeneout3_8h <- EBHMMNBMultiEM_2chain(Data=exprset3_8h, sizeFactors=sf, Conditions=conditions, UpdateRd=10, PriorFC=1.3)
GeneConfCalls3_8h <- GetConfidentCalls(EBHMMNBGeneout3_8h, FDR=.05, cutoff=.5, OnlyDynamic=T)
results3_8h <- as.data.frame(GeneConfCalls3_8h$Overall)
table(results3_8h$Most_Likely_Path)
EBHMMNBGeneout3_24h <- EBHMMNBMultiEM_2chain(Data=exprset3_24h, sizeFactors=sf, Conditions=conditions, UpdateRd=10, PriorFC=1.3)
GeneConfCalls3_24h <- GetConfidentCalls(EBHMMNBGeneout3_24h, FDR=.05, cutoff=.5, OnlyDynamic=T)
results3_24h <- as.data.frame(GeneConfCalls3_24h$Overall)
table(results3_24h$Most_Likely_Path)

## gene annotation ##
mart = useDataset("hsapiens_gene_ensembl", mart=useMart("ensembl"))
gSymbols_8h <- getBM(attributes =c('ensembl_gene_id','hgnc_symbol'), filters="ensembl_gene_id", values=rownames(results_8h), mart=mart)
gSymbols_24h <- getBM(attributes =c('ensembl_gene_id','hgnc_symbol'), filters="ensembl_gene_id", values=rownames(results_24h), mart=mart)
gSymbols2_8h <- getBM(attributes =c('ensembl_gene_id','hgnc_symbol'), filters="ensembl_gene_id", values=rownames(results2_8h), mart=mart)
gSymbols2_24h <- getBM(attributes =c('ensembl_gene_id','hgnc_symbol'), filters="ensembl_gene_id", values=rownames(results2_24h), mart=mart)
gSymbols3_8h <- getBM(attributes =c('ensembl_gene_id','hgnc_symbol'), filters="ensembl_gene_id", values=rownames(results3_8h), mart=mart)
gSymbols3_24h <- getBM(attributes =c('ensembl_gene_id','hgnc_symbol'), filters="ensembl_gene_id", values=rownames(results3_24h), mart=mart)

results_8h <- merge(results_8h, gSymbols_8h, by.x=0, by.y="ensembl_gene_id")
colnames(results_8h) <- c('ensembl_id','most_likely_path','max_PP','gene_symbol')
results_8h <- results_8h[!(results_8h$gene_symbol == ""), ]
write.csv(results_8h, "/media/shinigami/Elements/RNAseq/median_combat_ebseqhmm_8h.csv", row.names=F, quote=F)
results_24h <- merge(results_24h, gSymbols_24h, by.x=0, by.y="ensembl_gene_id")
colnames(results_24h) <- c('ensembl_id','most_likely_path','max_PP','gene_symbol')
results_24h <- results_24h[!(results_24h$gene_symbol == ""), ]
write.csv(results_24h, "/media/shinigami/Elements/RNAseq/median_combat_ebseqhmm_24h.csv", row.names=F, quote=F)

results2_8h <- merge(results2_8h, gSymbols2_8h, by.x=0, by.y="ensembl_gene_id")
colnames(results2_8h) <- c('ensembl_id','most_likely_path','max_PP','gene_symbol')
results2_8h <- results2_8h[!(results2_8h$gene_symbol == ""), ]
write.csv(results2_8h, "/media/shinigami/Elements/RNAseq/quantile_combat_ebseqhmm_8h.csv", row.names=F, quote=F)
results2_24h <- merge(results2_24h, gSymbols2_24h, by.x=0, by.y="ensembl_gene_id")
colnames(results2_24h) <- c('ensembl_id','most_likely_path','max_PP','gene_symbol')
results2_24h <- results2_24h[!(results2_24h$gene_symbol == ""), ]
write.csv(results2_24h, "/media/shinigami/Elements/RNAseq/quantile_combat_ebseqhmm_24h.csv", row.names=F, quote=F)

results3_8h <- merge(results3_8h, gSymbols3_8h, by.x=0, by.y="ensembl_gene_id")
colnames(results3_8h) <- c('ensembl_id','most_likely_path','max_PP','gene_symbol')
results3_8h <- results3_8h[!(results3_8h$gene_symbol == ""), ]
write.csv(results3_8h, "/media/shinigami/Elements/RNAseq/tmm_combat_ebseqhmm_8h.csv", row.names=F, quote=F)
results3_24h <- merge(results3_24h, gSymbols3_24h, by.x=0, by.y="ensembl_gene_id")
colnames(results3_24h) <- c('ensembl_id','most_likely_path','max_PP','gene_symbol')
results3_24h <- results3_24h[!(results3_24h$gene_symbol == ""), ]
write.csv(results3_24h, "/media/shinigami/Elements/RNAseq/tmm_combat_ebseqhmm_24h.csv", row.names=F, quote=F)
