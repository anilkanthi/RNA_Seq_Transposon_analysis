#DESeq2 - RNA-Seq with replicates

setwd("/Users/kantha01/Documents/Server/rna_seq/365/results")
library(DESeq2)
library(hexbin)
library(ggplot2)
library(RColorBrewer)


directory <- "/Users/kantha01/Documents/Server/rna_seq/365/results"
directory

sampleFiles <- grep("counts2",list.files(directory),value=TRUE)
sampleFiles

#sampleFiles<-sampleFiles[-c(4:6)]
#sampleFiles

------------------------------------------------------
sampleName <- sub("(*)_S.*","\\1",sampleFiles)
sampleCondition <- sub("(*)_S.*","\\1",sampleFiles)
sampleCondition[1]<-"365_D"
sampleCondition[2]<-"365_D"
sampleCondition[3]<-"365_D"
sampleCondition[4]<-"365_I"
sampleCondition[5]<-"365_I"
sampleCondition[6]<-"365_I"
sampleCondition[7]<-"365_ID"
sampleCondition[8]<-"365_ID"
sampleCondition[9]<-"365_ID"
sampleCondition[10]<-"365_UI"
sampleCondition[11]<-"365_UI"
sampleCondition[12]<-"365_UI"


sampleTable2 <- data.frame(sampleName = sampleFiles, fileName = sampleFiles,condition = sampleCondition)
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable2,directory = directory,design= ~ condition)

#matrix with data and all sampels
data<-as.data.frame(counts(dds))
data_element=data
repeats_all = read.table("repeats_all.txt")[,1]
repeats = as.character(read.table("repeats_unique_noRNA.txt")[,1])
data_repeats = NULL
data_repeats = t(sapply(repeats, function(re){
  repeats_rows = colSums(data[which(repeats_all == re),])
  data_repeats = rbind(data_repeats, repeats_rows)
}))
which(is.na(data_repeats))
rownames(data_repeats) = repeats
rowsums = rowSums(data_repeats)

samples = sampleName
#samples = c("365_D1","365_D2","365_D3", "365_ID1","365_ID2","365_ID3","365_UI1","365_UI2","365_UI3")
####
colData = data.frame(sampleCondition)
rownames(colData) = colnames(data_repeats)
colnames(colData) = c("condition")
colnames(data_repeats) = samples
condition = factor(sampleCondition)
#DESeq steps
dds = DESeqDataSetFromMatrix(countData=data_repeats,
                             colData = colData,
                             design = ~ condition)



dds = DESeq(dds)
norm<-counts(dds, normalize=TRUE)
write.csv(norm,file="normalized.365.PL.TE.txt",quote=FALSE)
rld <- rlog(dds)

plotPCA(rld, intgroup="condition")
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
write.table(assay(vsd),"365.VSD.TEs.txt",quote=FALSE)
library(ggplot2)

data <- plotPCA(rld, returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
pdf("365_by_condition_PCA.pdf",width=8,height=5,useDingbats = FALSE)
ggplot(data, aes(PC1, PC2, color=condition)) +
  geom_point(size=8) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(text = element_text(size=20))
dev.off()


pdf("365_by_replicate_PCA.pdf",width=8,height=5,useDingbats = FALSE)
ggplot(data, aes(PC1, PC2, color=name)) +
  geom_point(size=8) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(text = element_text(size=20))
dev.off()


#norm_counts = counts(dds, normalize=TRUE)
#norm_counts_avg = NULL
#norm_counts_avg = cbind(norm_counts_avg, rowMeans(norm_counts[,1:2]), rowMeans(norm_counts[,3:6]))
#colnames(norm_counts_avg) = c("365", "Ctrl")
#norm_counts_avg_log = log(norm_counts_avg+1, 10)
#comparisons
comparisons = rbind(c("365_D", "365_UI"),c("365_ID", "365_UI"))

for(i in 1:nrow(comparisons)){
  res = results(dds, c("condition", comparisons[i,1], comparisons[i,2]))
  res_sig = res[abs(res[,2]) > 1 & res[,6] < 0.05 & !is.na(res[,6]),]
  nrow(res_sig)
  assign(paste(comparisons[i,1], "_", comparisons[i,2], sep = ""), res)
  assign(paste(comparisons[i,1], "_", comparisons[i,2], "_sig_log2fc_1_FDR5", sep = ""), res_sig)
}

for(i in 1:nrow(comparisons)){
  res = results(dds, c("condition", comparisons[i,1], comparisons[i,2]))
  res_sig = res[abs(res[,2]) > 0.5 & res[,6] < 0.05 & !is.na(res[,6]),]
  nrow(res_sig)
  assign(paste(comparisons[i,1], "_", comparisons[i,2], sep = ""), res)
  assign(paste(comparisons[i,1], "_", comparisons[i,2], "_sig_log2fc_0.5_FDR5", sep = ""), res_sig)
}


write.csv(`365_D_365_UI`, file="365_D_365_UI_All_results.csv",quote=FALSE)
write.csv(`365_ID_365_UI`, file="365_ID_365_UI_All_results.csv",quote=FALSE)
write.csv(`365_ID_365_UI_sig_log2fc_1_FDR5`, file="365_ID_365_UI_sig_log2fc_1FDR5_25Changes.csv",quote=FALSE)
write.csv(`365_D_365_UI`, file="365_D_365_UI_All_results.csv",quote=FALSE)
