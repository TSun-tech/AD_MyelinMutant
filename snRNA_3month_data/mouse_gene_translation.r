#2021, Ting Sun
#in order to integrate mouse and human dataset
#gene symbols needs to be translated
#this is an example script of translating mouse genes to human genes using biomaRt
#dataset: myelin mutant, all cells

#for scRNA-seq objects
library(Seurat)
library(dplyr)

#for visualization|
library(cowplot)
library(ggplot2)
library(patchwork)

#for translation
library(biomaRt)

#organizingt sparse matrix
library(Matrix)

#read in source data
#Mouse<- readRDS("./2020Oct_MyelinMutants_all.rds")

#extract expression matrix and meta.data matrix
exp<-as.matrix(Mouse@assays$RNA@counts)
meta.data<-Mouse@meta.data

#extract all gene names
gene_trans<-as.data.frame(rownames(exp))

colnames(gene_trans)<-c("mouse_gene")

gene_trans$human_gene<-gene_trans$mouse_gene

#translate gene symbols which are avalible in biomaRt reference

musGenes <- rownames(exp)
# Basic function to convert mouse to human gene names
convertMouseGeneList <- function(x){
require("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , 
                 mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
humanx <- unique(genesV2[, 2])
# Print the first 6 genes found to the screen
print(head(humanx))
return(genesV2)
}

temp<-convertMouseGeneList(musGenes)

head(temp)
dim(temp)

#check length for translateble genes
length(unique(temp$MGI.symbol))
length(unique(temp$HGNC.symbol))

#map to original gene list in the translation column
gene_trans$human_gene<-as.character(gene_trans$human_gene)

for(i in 1:nrow(gene_trans)){
    position<-match(gene_trans$mouse_gene[i], temp$MGI.symbol)[1]
    gene_trans$human_gene[i] <- as.character(temp$HGNC.symbol[position])
}

#for non-translated genes, change gene format to human genes
gene_trans$human_gene[is.na(gene_trans$human_gene)]<-toupper(gene_trans$mouse_gene[is.na(gene_trans$human_gene)])

for ( k in 1:nrow(exp)){
    pos<-match(rownames(exp)[k],gene_trans$mouse_gene)
    rownames(exp)[k]<-gene_trans$human_gene[pos]
}

head(exp)
all(is.na(rownames(exp)) ==FALSE)

#recreate seurat object for further iintegration analysis
seu<-CreateSeuratObject(counts = exp, meta.data = meta.data)

seu








