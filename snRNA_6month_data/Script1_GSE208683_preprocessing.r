#GSE208683, analysis script 1
#data collected from snRNA-seq of mouse brain hemispheres
#in total 4 genotypes: WT, CnpKO, 5xFAD, CnpKO 5xFAD

#this script contains prepsocessing of CellRanger aligned, CellBender processed UMI count matrices
#until unbiased clustering, cluster clean and cell type annotations

#run under sc_update, R 4.1.2, Seurat 4.1.1

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

#read in CellBender resuls (.h5)
dir<-list.dirs(path = "/.../RAW/CellBender_outs", 
               recursive = FALSE)

#prepare genotype details (name of samples) according to the sample list
Genotype<-c(rep("WT",2), rep("CnpKO",2), rep("AD", 2), rep("CnpKOAD", 2))
Sample_ID<-c("WT_1","WT_2","CnpKO_1","CnpKO_2","AD_1","AD_2","CnpKOAD_1","CnpKOAD_2")

(orig.id<-paste0(c(1:8),"_HS"))

#prepare list of matrix and read in data using Read10x
number<-8
expression_mat<-vector(mode = 'list', length = number)
names(expression_mat)<-orig.id

for (i in 1:8){
    expression_mat[[i]]<-Read10X_h5(filename = paste0(dir[[i]],"/",orig.id[i],"_output_filtered.h5"))
}

#organize seurat object list
seurat_list<-vector(mode = 'list', length = number)

#write in sample meta.data into seurat onjects
for (k in 1:8){
    seurat_list[[k]]<-CreateSeuratObject(counts = expression_mat[[k]], project = orig.id[[k]])
    seurat_list[[k]]@meta.data$Genotype<-as.character(Genotype[[k]])
    seurat_list[[k]]@meta.data$SampleID<-as.character(Sample_ID[[k]])
}

head(seurat_list[[1]]@meta.data)

#plot individual sample sequencing depth and define filter paramters for each sample
for (i in 1:8){
    seurat_list[[i]][["percent.mt"]]<-PercentageFeatureSet(seurat_list[[i]],pattern = "mt-")
    print(VlnPlot(seurat_list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                  ncol = 3, pt.size = 0.0001))
    print(FeatureScatter(seurat_list[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
}

#based on each plot, a filter value (high end) for individual sample is generated as an external csv
#column names for the data.fram are: orig.id, nFeature_RNA, nCount_RNA

#read in filter levels
filter<-read.csv("...indir/Seurat_raw_data_individual_sample_filter.csv")

#before filter, observe defined cutoffs on the plot
for(m in 1:8){
    print(FeatureScatter(seurat_list[[m]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+
    geom_hline(yintercept = c(300,filter[m,2]), linetype = "dashed")+
    geom_vline(xintercept = c(500,filter[m,3]), linetype = "dashed"))
    
    #back up plots
    pdf(file = paste0(".../outdir/",
                     orig.id[[m]],"_cell_filter.pdf"))
    print(FeatureScatter(seurat_list[[m]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+
    geom_hline(yintercept = c(300,filter[m,2]), linetype = "dashed")+
    geom_vline(xintercept = c(500,filter[m,3]), linetype = "dashed"))
    dev.off()
}


#after confirmation the cutoffs, apply for individual sample filter
#replot the post-filtered data
for (k in 1:8){
    seurat_list[[k]]<-subset(seurat_list[[k]], subset = nFeature_RNA > 300 & nFeature_RNA < filter[k,2] & nCount_RNA > 500 & nCount_RNA < filter[k,3])
    
    print(FeatureScatter(seurat_list[[k]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+
    geom_hline(yintercept = c(300,filter[k,2]), linetype = "dashed")+
    geom_vline(xintercept = c(500,filter[k,3]), linetype = "dashed"))
}

#merge object list into one seurat object

seu<-merge(seurat_list[[1]],y = c(seurat_list[2:8]), project = "DeppSun_2022", add.cell.ids = Sample_ID)

seu

#observe if merged object has obvious outlier cells
print(FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))

#proceed with data normalization, scaling, and linear dimentionality reduction

seu<-NormalizeData(seu)

seu<-FindVariableFeatures(seu, selection.method = 'vst',
                            nfeatures = 3000)

all.genes <- rownames(seu)
seu <- ScaleData(seu, features = all.genes)

seu<-RunPCA(seu, features = VariableFeatures(object = seu), npcs = 100)
ElbowPlot(seu, ndims = 100)

#perform neighbouring embedding using UMAP
seu<-RunUMAP(seu, dims = 1:30)

DimPlot(seu, group.by = "SampleID", reduction = "umap")
DimPlot(seu, group.by = "Genotype", reduction = "umap")

#confirmed no batch effects between samples
#perform clustering analysis and calculate cluster marker gene
#for further clean up of doublets and annotating cell populations
set.seed(20220215)

seu <- FindNeighbors(object = seu, dims = 1:30)

seu <- FindClusters(object = seu, resolution = 0.8)
table(Idents(seu))

DimPlot(object = seu, reduction = 'umap',label = TRUE, 
        pt.size = 0.1,label.size = 6, repel = T) + theme(aspect.ratio=1)+NoLegend()

#apply default marker gene cauculation, threshold min.pct =0.25, logfc 0.25
#calculate marker genes using RNA assay
DefaultAssay(seu)<-"RNA"
Idents(seu)<-"RNA_snn_res.0.8"

markers_rna<-FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#export marker genes for annotating clusters
#based on marker gene calculations, cells are allocated as specific celltypes, or noise populations/doublet
#cleaned object expression matrix and corresponding meta.data are uploaded to GSE208683








































