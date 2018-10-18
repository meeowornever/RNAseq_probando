
library(EBSeqHMM)
library(EBSeq)
library(sva)
library(edgeR)
library(biomaRt)
library(RColorBrewer)

## Loading dataset
aitordata <- read.table("/media/shinigami/Elements/RNAseq/countsaitor.txt", comment.char="#")
aitordata <- as.matrix(aitordata)
rownames(aitordata) <- aitordata[,1]
aitordata <- aitordata[,-1]
aitordata <- aitordata[-1,]
exprMatrix <- aitordata[,6:29]
colnames(exprMatrix) <- c("D8R1", "D24R3", "D24R2", "D8R3", "D8R2", "D24R1", "M24R3", "M8R2", "M8R1", "M8R3", "M24R2", "M24R1", "T8R1", "T24R2", "T24R3", "T24R1", "T8R2", "T8R3", "TM8R3", "TM8R1", "TM24R2", "TM24R1", "TM8R2", "TM24R3")
mode(exprMatrix) <- 'numeric'

## Filtering of low-abundant genes
keep <- filterByExpr(exprMatrix)
exprMatrixflt <- exprMatrix[keep,]
checker <- aveLogCPM(exprMatrixflt)
expressed <- exprMatrixflt[checker>1,]

##phenotype matrix
samples_series <- colnames(expressed)
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

## Normalization (Median)
sizes <- MedianNorm(expressed)
GeneNormData <- GetNormalizedMat(expressed, sizes)

boxplot(log2(GeneNormData), las=2, col=col, main="")
title(main="MedianNorm Normalised data",ylab="Log2-counts")
plotMDS(log2(GeneNormData), top = 1000, labels = NULL, col = as.numeric(phenodata$condition), cex = 2)
title(main="Treatment groups")
plotMDS(log2(GeneNormData), top = 1000, labels = NULL, col = as.numeric(phenodata$replicate), cex = 2)
title(main="Replicate groups")

## Batch Effect Correction (Combat)
logCPM <- cpm(GeneNormData, log=T, prior.count=.5)
batch <- phenodata$replicate
modcombat <- model.matrix(~condition, data=phenodata)
cleaned <- ComBat(dat=logCPM, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

boxplot(cleaned, las=2, col=col, main="")
title(main="MedianNorm+Combat data",ylab="Log2-counts")
plotMDS(cleaned, top = 1000, labels = NULL, col = as.numeric(phenodata$condition), cex = 2)
title(main="Treatment groups")
plotMDS(cleaned, top = 1000, labels = NULL, col = as.numeric(phenodata$replicate), cex = 2)
title(main="Replicate groups")

antilog<-function(lx,base) {
  lbx<-lx/log(exp(1),base=base)
  result<-exp(lbx)
  result
}
cleanedanti <- antilog(cleaned, 2)
exprset_8h <- cleanedanti[,phenodata$timepoint == 8]
exprset_24h <- cleanedanti[,phenodata$timepoint == 24]

## EBSEQHMM
sf <- rep(1,12)
conditions <- as.factor(rep(c('dmso','mkc','taxol','taxolmkc'), each=3))

EBHMMNBGeneout8h_sva <- EBHMMNBMultiEM_2chain(Data=exprset_8h, sizeFactors=sf, Conditions=conditions, UpdateRd=10, PriorFC=1.3)
GeneDECalls <- GetDECalls(EBHMMNBGeneout8h_sva, FDR=.05)
GeneConfCalls8h_sva <- GetConfidentCalls(EBHMMNBGeneout8h_sva, FDR=.05, cutoff=.5, OnlyDynamic=T)
results_8h <- as.data.frame(GeneConfCalls8h_sva$Overall)
table(results_8h$Most_Likely_Path)

EBHMMNBGeneout24h_sva <- EBHMMNBMultiEM_2chain(Data=exprset_24h, sizeFactors=sf, Conditions=conditions, UpdateRd=10, PriorFC=1.3)
GeneConfCalls24h_sva <- GetConfidentCalls(EBHMMNBGeneout24h_sva, FDR=.05, cutoff=.48, OnlyDynamic=T)
results_24h <- as.data.frame(GeneConfCalls24h_sva$Overall)
table(results_24h$Most_Likely_Path)

## Gene annotation
mart = useDataset("hsapiens_gene_ensembl", mart=useMart("ensembl"))
gSymbols8h <- getBM(attributes =c('ensembl_gene_id','hgnc_symbol'), filters="ensembl_gene_id", values=rownames(results_8h), mart=mart)
gSymbols24h <- getBM(attributes =c('ensembl_gene_id','hgnc_symbol'), filters="ensembl_gene_id", values=rownames(results_24h), mart=mart)

results_8h <- merge(results_8h, gSymbols8h, by.x=0, by.y="ensembl_gene_id")
colnames(results_8h) <- c('ensembl_id','most_likely_path','max_PP','gene_symbol')
results_8h <- results_8h[!(results_8h$gene_symbol == ""), ]
write.csv(results_8h, "/media/shinigami/Elements/RNAseq/res8h_combat.csv", row.names=F, quote=F)

results_24h <- merge(results_24h, gSymbols24h, by.x=0, by.y="ensembl_gene_id")
colnames(results_24h) <- c('ensembl_id','most_likely_path','max_PP','gene_symbol')
results_24h <- results_24h[!(results_24h$gene_symbol == ""), ]
write.csv(results_24h, "/media/shinigami/Elements/RNAseq/res24h_combat.csv", row.names=F, quote=F)
