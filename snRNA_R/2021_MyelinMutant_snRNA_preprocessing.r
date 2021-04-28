#2020Dec, Ting Sun
#myelin mutants sample 1-8 (all samples)
#organizationg of original data into Seurat object
#script include read10X matrix
#form individual seurat object and apply filter
#save individual object
#merge samples into one temporary object, run until clustering for preliminary analysis

#run under sc_env


#attach packages

#for scRNA-seq objects
library(Seurat)
library(dplyr)
library(Matrix)
library(abind)
library(sctransform)

#for parallelization
library(future)
library(future.apply)

#for visualization
library(cowplot)
library(ggplot2)

#for DGE analysis
library(DESeq2)
library(MAST)
library(patchwork)
library(pheatmap)

#trajectory analysis
library(scater)
library(slingshot)
library(ggbeeswarm)

#support general object transfer
library(SingleCellExperiment)

#normally unnecessary, for loading sparse matrix
#library(DropSeq.util)

#read in CellRanger aligned results
data_dir<-'/media/tsun/Data/Tsun/Ting/Myelin_mutant/Aligned/filtered_feature/'

list.files(data_dir)
folders<-list.files(data_dir)
length(folders)

#label sample genotypes and replicate ID
genotype<-c("WT","CnpKO","Plp1KO","Foxg1MbpKO","WT","CnpKO","Plp1KO","Foxg1MbpKO")
sample_ID<-c("WT1","CnpKO1","Plp1KO1","Foxg1MbpKO1",
            "WT2","CnpKO2","Plp1KO2","Foxg1MbpKO2")

#prepare empty lists to store constructed seurat objects
number<-8
expression_mat<-vector(mode = 'list', length = number)
seurat_list<-vector(mode = 'list', length = number)


for (i in 1:8){
    temp<-paste0('/media/tsun/Data/Tsun/Ting/Myelin_mutant/Aligned/filtered_feature/',
                                                folders[[i]])
    list.files(temp)
    expression_mat[[i]]<-Read10X(data.dir = temp)
}

for (k in 1:8){
    seurat_list[[k]]<-CreateSeuratObject(counts = expression_mat[[k]], project = folders[[k]])
    seurat_list[[k]]@meta.data$Genotype<-as.character(genotype[[k]])
    seurat_list[[k]]@meta.data$SampleID<-as.character(sample_ID[[k]])
}

seurat_list

#filter individual dataset based on sequencing depth and aligned gene numbers
for (i in 1:8){
    print(VlnPlot(seurat_list[[i]], features = c("nFeature_RNA", "nCount_RNA"), ncol = 2))
    print(FeatureScatter(seurat_list[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
}

seurat_list[[1]]<-subset(seurat_list[[1]],subset = nFeature_RNA > 200 & nFeature_RNA < 9400 & nCount_RNA > 500 & nCount_RNA < 75000)
seurat_list[[2]]<-subset(seurat_list[[2]],subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & nCount_RNA > 500 & nCount_RNA < 75000)
seurat_list[[3]]<-subset(seurat_list[[3]],subset = nFeature_RNA > 200 & nFeature_RNA < 9500 & nCount_RNA > 500 & nCount_RNA < 90000)
seurat_list[[4]]<-subset(seurat_list[[4]],subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & nCount_RNA > 500 & nCount_RNA < 100000)
seurat_list[[5]]<-subset(seurat_list[[5]],subset = nFeature_RNA > 200 & nFeature_RNA < 11000 & nCount_RNA > 500 & nCount_RNA < 100000)
seurat_list[[6]]<-subset(seurat_list[[6]],subset = nFeature_RNA > 200 & nFeature_RNA < 11000 & nCount_RNA > 500 & nCount_RNA < 130000)
seurat_list[[7]]<-subset(seurat_list[[7]],subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & nCount_RNA > 500 & nCount_RNA < 130000)
seurat_list[[8]]<-subset(seurat_list[[8]],subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & nCount_RNA > 500 & nCount_RNA < 100000)

#attach sample ID before cell barcode

sample_ID

seu<-merge(seurat_list[[1]],y = c(seurat_list[2:8]), project = "MyelinMutant", add.cell.ids = sample_ID)

seu

#head(seu@meta.data)

#perform normalization and scaling
seu<-NormalizeData(seu)

seu<-FindVariableFeatures(seu, selection.method = 'vst',
                            nfeatures = 3000)

all.genes <- rownames(seu)
seu <- ScaleData(seu, features = all.genes)

#saveRDS(seu,'/media/tsun/Data/Tsun/Ting/Myelin_mutant/RDS/2020Oct_MyelinMutants_all.rds')










































