#2020, Ting Sun
#for Depp and Sun et al manuscript
#microglia bulk RNA sequencing
#raw fastq files aligned to reference genome using STAR default parameters
#gene counts extracted by featureCounts

#sequenced genotypes (6-month-old)
#WT, CnpKO, 5xFAD, CnpKO x 5xFAD

#script aim to summarize raw counts and TPM counts table
#use normalized value for PCA embedding and verify potential outlier
#extract PC1 most relevant genes to obtain global view of gene regulation between all genotypes


library(DESeq2)
library(dplyr)
library(ggplot2)

#for visualization
    library(ggplot2)
    library(geneplotter)
    library(RColorBrewer)
    #library(pheatmap)

    library(dplyr)
    library(tidyr)

sessionInfo()

#first read in aligned gene counts files and summary raw counts table
count.list<-list.files(path = "./",
                       pattern = "tabular")

count.list

meta.data<-matrix(ncol = 3, nrow = 16)

colnames(meta.data)<-c("file_name","seq_number","genotype")

meta.data<-as.data.frame(meta.data)

meta.data$file_name<-count.list

for(i in 1:16){
    temp<-strsplit(meta.data$file_name[i],"833s")[[1]][2]
    meta.data$seq_number[i]<-as.numeric(strsplit(temp,"_")[[1]][1])
}

meta.data

meta.data$seq_number<-as.numeric(meta.data$seq_number)

meta.data<-meta.data[order(meta.data$seq_number, decreasing = F),]

meta.data

meta.data$genotype<-c(rep("WT",4),rep("Cnp_ko",4),rep("5xFAD",4),rep("Cnp_5xFAD",4))

meta.data$sample<-NA

meta.data$sample<-paste0(meta.data$genotype,"_",meta.data$seq_number)

meta.data

#combine gene raw counts
readtab<-function(x){
    temp<-read.csv(x, sep = "",stringsAsFactors = FALSE)
    rownames(temp)<-temp$Geneid
    temp<-temp[-1]
    return(temp)
}

setwd("./")
counts.files<-lapply(count.list, readtab)

str(counts.files)

raw<-Reduce(merge, lapply(counts.files, function(x) data.frame(x, rn = row.names(x))))

head(raw)

rownames(raw)<-raw$rn
raw<-raw[,-1]

colnames(raw)

for (i in 1:16){
    position<-grep(colnames(raw)[[i]], meta.data$file_name)
    colnames(raw)[[i]]<-meta.data$sample[position]
    
}

colnames(raw)

raw<-raw[,c(1,9:16,2:4,5:8)]

head(raw)
colnames(raw)

#translate gene ID to gene symbol
library("org.Mm.eg.db")

essemble<-rownames(raw)

    convert<-mapIds(org.Mm.eg.db, keys = essemble, keytype = "ENSEMBL", column="SYMBOL")

head(convert)
length(convert)

library(biomaRt)
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
annot<-getBM(c("ensembl_gene_id", "mgi_symbol", "chromosome_name", "strand", "start_position", "end_position","gene_biotype"), 
             mart=ensembl)

head(annot)
dim(annot)

raw$gene_symbol<-NA

for(i in 1:nrow(raw)){
    position<-match(as.character(rownames(raw)[[i]]),annot$ensembl_gene_id)
    raw$gene_symbol[[i]]<-annot$mgi_symbol[position]
}

raw<-raw[,c(17,1:16)]

head(raw)



#caulate TPM value

#read in gene length matrix

len<-read.csv("/media/tsun/Data/Tsun/Constanze/2020_microglia_RNAseq/feature_length/featureCounts on collection 58: Feature lengths/p833s1_Nave_S42_L005_R1_001.fastq.gz.tabular",
             sep="\t",stringsAsFactors=FALSE)


## functions for tpm
## from https://gist.github.com/slowkow/c6ab0348747f86e2748b#file-counts_to_tpm-r-L44
## from https://www.biostars.org/p/171766/

tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

#prepare empty matrix
res<-matrix(ncol = 19, nrow = nrow(test))

colnames(res)[1:3]<-c("ENS_ID","Gene_symbol","Gene_length")
colnames(res)[4:19]<-c(paste0(colnames(raw)[2:17],"_TPM"))

res$ENS_ID<-as.character(raw$X)
res$Gene_symbol<-raw$gene_symbol

for(k in 1:nrow(res)){
    position<-match(res$ENS_ID[k],len$Geneid)
    res$Gene_length[k]<-len$Length[position]
}

for (m in 4:19){
    res[,m]<-tpm(raw[,(m-1)],res$Gene_length)
}
#res == TPM value matrix



#use tpm values for PCA embedding to detect if there are outliers
tpm_pca<-res
rownames(tpm_pca)<-tpm_pca$ENS_ID
tpm_pca<-tpm_pca[,-c(1:4)]
head(tpm_pca)

tpm_pca_cal<-prcomp(t(tpm_pca), scale = F)
str(tpm_pca_cal)

#extract sd and PC1, PC2 information, construct new matrix for visualization
percentVar <- round(100*tpm_pca_cal$sdev^2/sum(tpm_pca_cal$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = tpm_pca_cal$x[,1], PC2 = tpm_pca_cal$x[,2],
                    genotype = metadata$genotype)

#plot PCA
(manuscript_plot<-ggplot(dataGG, aes(PC1, PC2)) +
      geom_point(aes(shape = genotype, colour = genotype)) +
  ggtitle("PCA plot of Microglia RNA-seq TPM profiles") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio) +
#first change shape to dots
  scale_shape_manual(values = c(rep(16,4))) + 
#color change based on RGB values
  scale_color_manual(values = c("#76069AFF",'#F8B914FF', '#48ADF0FF',
                               "gray33")))

#extract top 100 genes of PC1
var<-get_pca_var(pca)

pc1<-var$contrib[,c("Dim.1")]

pc1<-sort(pc1, decreasing = T)
pc1<-as.data.frame(pc1)

#use gene annotation to translate gene ID to gene symbol








