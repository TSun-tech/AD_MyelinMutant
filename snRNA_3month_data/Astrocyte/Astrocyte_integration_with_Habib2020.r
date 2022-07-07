#2021, Ting Sun
#myelin mutant snRNA-seq
#Astrocyte subset
#only CnpKO, MbpKO and WT
#integration analysis with DAA (Habib 2020)

#run under sc_env

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

ast<-readRDS("./2021March_MyelinMutants_Allsamples_Astrocyte.rds")
daa<-readRDS("./GSE143758_Admouse_Hippocampus_7m_AllNuclei.rds")

ast
daa

#organize objects to have matched meta.data

#subset objects
#myelin mutant: CnpKO, Foxg1MbpKO, WT
#DAA dataset

ast<-subset(ast, subset = Genotype %in% c("CnpKO","Foxg1MbpKO","WT"))

all<-merge(ast, daa, project = "Depp_Sun2021")

all

head(all@meta.data)

#primary embedding analysis
#detect potential batch effect
#one step normalisation until PCA embedding

all<-NormalizeData(all)

all<-FindVariableFeatures(all, selection.method = 'vst',
                            nfeatures = 2000)

all.genes <- rownames(all)
all <- ScaleData(all, features = all.genes)

all<-RunPCA(all, features = VariableFeatures(object = all))
ElbowPlot(all, ndims = 50)

all<-RunUMAP(all, dims = 1:20)

DimPlot(all, group.by = "Study")
DimPlot(all, group.by = "Genotype")
DimPlot(all, group.by = "CellType")
DimPlot(all, group.by = "SampleID", label = T)+NoLegend()

#identify study as batch effect
#also correct for both datasets sequencing batches to obtain better integration

#seperate batch 1 and batch2 from myelin mutant data before processing with integration pipeline
all@meta.data$batch<-all@meta.data$Study

all@meta.data$batch[which(all@meta.data$Study == "MyelinMutant")]<-all@meta.data$SampleID[which(all@meta.data$Study == "MyelinMutant")]
all@meta.data$batch[which(all@meta.data$Study == "Habib2020")]<-all@meta.data$treatment[which(all@meta.data$Study == "Habib2020")]


unique(all@meta.data$batch)

all@meta.data$batch<-factor(all@meta.data$batch,
                           levels = unique(all@meta.data$batch),
                           labels = c(rep("batch1",3), rep("batch2",3),
                                      "EZ","NP40"))
unique(all@meta.data$batch)

#prepare to correct for sample batch effect
#SCTranform pipeline

DefaultAssay(all)<-"RNA"
seu.list<-SplitObject(all, split.by = "batch")

for (k in c(1:length(seu.list))) {
    seu.list[[k]] <- SCTransform(seu.list[[k]], verbose = TRUE)
}

seu.features <- SelectIntegrationFeatures(object.list = seu.list, nfeatures = 1000)
seu.list <- PrepSCTIntegration(object.list = seu.list, anchor.features = seu.features)

seu.anchors <- FindIntegrationAnchors(object.list = seu.list, normalization.method = "SCT", 
                                        anchor.features = seu.features
                                      #, reference = c(3) #choose Habib2020 as reference dataset
                                     )

seu.integrated <- IntegrateData(anchorset = seu.anchors, normalization.method = "SCT")

seu.integrated <- RunPCA(seu.integrated, verbose = FALSE, npcs = 50)
ElbowPlot(seu.integrated, ndims = 50)

seu.integrated <- RunUMAP(seu.integrated, dims = 1:30)

#evaluate integration result

DimPlot(seu.integrated, group.by = "Study")
DimPlot(seu.integrated, group.by = "Genotype")
DimPlot(seu.integrated, group.by = "orig.cluster", label = TRUE)
DimPlot(seu.integrated, group.by = "Age")
DimPlot(seu.integrated, group.by = "SampleID", label = TRUE)+NoLegend()

#perform neibouring embedding to detect potential DAA population
seu.integrated<-FindNeighbors(seu.integrated, dims = 1:30)

seu.integrated <- FindClusters(object = seu.integrated, resolution = 0.5)
table(Idents(seu.integrated))

DimPlot(object = seu.integrated, reduction = 'umap',label = TRUE, 
        pt.size = 0.1,label.size = 6, repel = T) + theme(aspect.ratio=1)

#highlight each genotype to check integration quality and potential DAA cluster
seu.integrated@meta.data$Genotype[which(seu.integrated@meta.data$Genotype == "WT")] <- paste0(seu.integrated@meta.data$Study[which(seu.integrated@meta.data$Genotype == "WT")],
                                                                                             "_","WT")

unique(seu.integrated@meta.data$Genotype)


#original genotype highlight
color<- c("orange","#749dae","#E07882",'plum4',"tan")

sample_names<-unique(seu.integrated@meta.data$Genotype)

sampleID_list<-vector(mode = "list", length = length(sample_names))

Idents(seu.integrated)<-"Genotype"

for(i in 1:length(sample_names)){
    sampleID_list[[i]]<-WhichCells(seu.integrated, idents = sample_names[i])
}

for (k in 1:length(sample_names)){
        print(DimPlot(seu.integrated, label=F, group.by="Genotype", 
        cells.highlight= sampleID_list[[k]], pt.size = 0.01)+ 
              scale_color_manual(labels = c("Other Genotype", as.character(sample_names[k])), 
                                 values = c("grey", color[k])) +
              labs(color = "legend title"))
    }

(plt1<-DimPlot(seu.integrated, group.by = "Genotype",
              cols = c('plum4',"#749dae","#E07882","tan","orange")))

#hilighting plots
#for Habib cluster4 (DAA) on integrated object, red
Idents(seu.integrated)<-"orig.cluster"
(daa_umap<-DimPlot(seu.integrated, label=F, group.by="orig.cluster", 
        cells.highlight= WhichCells(seu.integrated, idents = "Habib2020_4"))+ 
              scale_color_manual(labels = c("Other cluster", "DAA"), 
                                 values = c("grey", "lightseagreen")) +
              labs(title = "Highlight UMAP for DAA cluster"))












