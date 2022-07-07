#2012April06, Ting Sun
#myelin mutant scnRNA-seq
#figure preperation for manuscript
#microglia and OLG check point genes over different cell subpopulations
#use same y-axis to make sure comparable plot

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

#ABA spatial map
library(voxhunt)

#read in datasets
#mic<-readRDS("./2020Oct_MyelinMutants_Allsamples_Microglia.rds")
#olg_all<-readRDS("./2020Oct_MyelinMutants_Allsamples_Oligodendrocyte_run2.rds")

mic
olg_all

#first subset and recalculate objects to only CnpKO and WT
mic<-subset(mic, subset = Genotype %in% c("CnpKO", "WT"))
olg<-subset(olg_all, subset = Genotype %in% c("CnpKO", "WT"))

table(olg@meta.data$Genotype)

#combine two cell types
all<-merge(olg, mic)

all

#one step normalisation until PCA embedding

all<-NormalizeData(all)

all<-FindVariableFeatures(all, selection.method = 'vst',
                            nfeatures = 2000)

all.genes <- rownames(all)
all <- ScaleData(all, features = all.genes)

all<-RunPCA(all, features = VariableFeatures(object = all))
ElbowPlot(all, ndims = 50)

all<-RunUMAP(all, dims = 1:15)

DimPlot(all, reduction = "umap", group.by = "Genotype")
DimPlot(all, reduction = "umap", group.by = "CellType")

markers<-c("Pdgfra","Cspg4",#OPC markers,
          "Enpp6","Bcas1", #NFO markers
          "Mbp","Plp1","Mog","Mag",#MOL markers
          "Cx3cr1","P2ry12","Tmem119",
        "Trem2","Apoe","Apod","Tyrobp","Ms4a7","Fth1",
          "Acsl1","Dpyd")

#further clean up
all<-subset(all, subset = CellType == "unknown", invert = T)
all<-NormalizeData(all)

unique(all@meta.data$CellType)

all@meta.data$CellType<-factor(all@meta.data$CellType,
                              levels = c("OPC","Newly formed Oligodendrocyte",
                                        "Mature Oligodendrocyte",
                                        "Homeostatic","Myelin processing",
                                        "DAM", "Brain Border \n Marcrophages"),
                              labels = c("OPC","Newly Formed Oligodendrocyte",
                                        "Mature Oligodendrocyte",
                                        "Homeostatic Microglia",
                                         "Myelin Processing Microglia",
                                        "DAM", "Brain Border Marcrophages"))

Idents(all)<-"CellType"
(plt<-VlnPlot(all, features = markers, stack = TRUE, 
        group.by = "CellType",
       flip = T,
       fill.by = "ident" #, slot = "scale.data"
       ,cols = c("bisque3","darksalmon","darkolivegreen3",
                 "tan1","lightskyblue","lightcoral","lightseagreen")
       )+NoLegend())






















