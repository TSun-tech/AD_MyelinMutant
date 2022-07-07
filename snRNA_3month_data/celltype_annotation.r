#2020, Ting Sun
#myelin mutants all samples
#preliminary analysis
#clustering and cell type annotation
#cell quantity and quality check
#cell proportion check
#first analysis of microglia Trem2 expression potential differences among myelin mutants

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

#read in preprocessed Seurat object
seu<-readRDS("./2020Oct_MyelinMutants_all.rds")

#perform dimensionality reduction using PCA with 3000 calculated variable features

seu<-RunPCA(seu, features = VariableFeatures(object = seu))

#check suitable number of PCs for neighbouring embedding
ElbowPlot(seu, ndims = 50)

#use first 50 PCs for first round neighbouring embedding (UMAP) analysis
seu<-RunUMAP(seu, dims = 1:50)

#perform unbiased clustering analysis
#find neibours
#set.seed(20201214)

seu <- FindNeighbors(object = seu, dims = 1:50)

#after tests, use resolution 0.5 for determine clusters

seu <- FindClusters(object = seu, resolution = 0.5)
table(Idents(seu))

#visualize cluster location
DimPlot(object = seu, reduction = 'umap',label = TRUE, 
        pt.size = 0.1,label.size = 6, repel = T) + theme(aspect.ratio=1)+NoLegend()

#calculate cluster marker genes
#threshold pct > 0.25, LogFC > 0.25
DefaultAssay(seu)<-"RNA"
markers_rna<-FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

write.csv(markers_rna,"./2020_MyelingMutants_scRNA_allSamples_res05_cluster_marker_RNA.csv")

#chek marker gene table of each cluster in detail and annotate cluster cell type identity
#mark and remove noisy clusters (potential coublet or apoptotic cells)

#before data clean up, visualize CNS cell type marker genes

GABA<-c("Rbfox3","Gad1","Gad2","Gadd45b")
Gluta<-c("Grin2b","Meis2","Gls")
Gly<-c("Slc6a9")
#NPC<-("Nes","Sox1","Smarca4","Pax3","Pax6")
MOL<-c("Mbp","Cnp","Mog","Olig1","Olig2","Plp1","Car2")
NFOL<-c("Sirt2","Bcas1","Bcan")
OPC<-c("Cspg4","Pdgfra","Sox6","Itpr2")
astro<-c("Gfap","S100b","Aqp4")
microglia<-c("Cx3cr1","P2ry12","Cd68","Itgam","Trem2","Hexb","Aif1")

cell_gene<-c(GABA,Gluta,Gly,MOL,NFOL,OPC,astro,microglia)
cell_gene

for( k in cell_gene){
    print(FeaturePlot(seu, features = k))
}

#call cluster identity, combine with marker gene calculation
#and typical marker visualization
#define and annotate major cell types and filter outlier cells

Idents(seu)<-"seurat_clusters"

DimPlot(object = seu, reduction = 'umap',label = TRUE, 
        pt.size = 0.1,label.size = 6, repel = T) + theme(aspect.ratio=1)+NoLegend()

seu@meta.data$CellType<-seu@meta.data$seurat_clusters

seu@meta.data$CellType<-factor(seu@meta.data$CellType,
                              levels = c(0:35),
                              labels = c("Mature Oligodendrocyte","Neuron","Neuron","Astrocyte","Neuron","Neuron",
                                        "unknown","Neuron","Inhibitory neuron","Microglia","OPC",
                                        "Inhibitory neuron","Endothelial","Neuron","Neuron","Neuron",
                                        "unknown","Inhibitory neuron","Neuron","Neuron","Neuron",
                                        rep("Neuron",3),"Endothelial","Pericyte",
                                        "Neuron","unknown","unknown","Astrocyte","Neuron",
                                        "Newly formed Oligodendrocyte", "Microglia","unknown","Neuron","unknown"))

DimPlot(seu,reduction = "umap",group.by = "CellType", label = TRUE)+NoLegend()

seu<-subset(seu,subset = CellType == "unknown", invert = TRUE)

seu

DimPlot(seu,reduction = "umap",group.by = "CellType", label = TRUE)+NoLegend()

#reapply normalization, scaling, dimentionality reduction and neighbouring embedding
#to cleaned up object
seu<-NormalizeData(seu)

seu<-FindVariableFeatures(seu, selection.method = 'vst',
                            nfeatures = 3000)

all.genes <- rownames(seu)
seu <- ScaleData(seu, features = all.genes)

seu<-RunPCA(seu, features = VariableFeatures(object = seu))
ElbowPlot(seu, ndims = 50)

seu<-RunUMAP(seu, dims = 1:50)

DimPlot(seu, group.by = "SampleID")
DimPlot(seu, group.by = "CellType")

#visualize marker genes after re-apply embedding to cleaned object
for( k in cell_gene){
    print(FeaturePlot(seu, features = k))
}
























