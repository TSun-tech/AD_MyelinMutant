#2021, Ting Sun
#Myelin mutant scRAN-seq data
#Microglia subset

#merging with external datasets

#in house: Myelin mutant microglia profile
#Genotype: CnpKO with corresponding WT

#External dataset
#AD: Zhou 2020 7-month (5xFAD)

#chcek underlying batch effect and perform SCTransform pipeline for elimination
#embedding and clustering analysis using integrated assay
#goal of analysis: annotate DAM and myelin processing population among myelin mutant and AD mice
#use Pseudotime trajectory analysis to dissect microglia activation stage differences in CnpKO vs other genotypes

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

#read in source data
#in house myelin mutant OLG subset
##ting<-readRDS("./2020Oct_MyelinMutants_Allsamples_Microglia.rds")

#Zhou et.al 2020 (5xFAD 7month)
#zhou<-readRDS("./mouse_7m_samples_sclaed_celltype_identified.rds")

ting
zhou

#subset microglia from external dataset
zhou<-subset(zhou, subset = CellType == "microglia")

#check all meta.data colnames after organization
colnames(ting@meta.data)
colnames(zhou@meta.data)



#dataset merging
obj<-list(ting, zhou)
names(obj)<-c("Myelin Mutant", "5xFAD")

mic<-merge(x = obj[[1]], y = c(obj[2]), project = "Depp_Sun2021")

#focus on CnpKO with corresponding WT and AD genotypes
mic<-subset(mic, subset = Genotype %in% c("WT","CnpKO","WT_5xFAD","Trem2KO","Trem2KO_5xFAD"))

#primary embedding analysis
#detect potential batch effect
#one step normalisation until PCA embedding

mic<-NormalizeData(mic)

mic<-FindVariableFeatures(mic, selection.method = 'vst',
                            nfeatures = 2000)

all.genes <- rownames(mic)
mic <- ScaleData(mic, features = all.genes)

mic<-RunPCA(mic, features = VariableFeatures(object = mic))
ElbowPlot(mic, ndims = 50)

mic<-RunUMAP(mic, dims = 1:20)

DimPlot(mic, group.by = "Study")
DimPlot(mic, group.by = "Genotype")
DimPlot(mic, group.by = "CellType")
DimPlot(mic, group.by = "SampleID", label = T)+NoLegend()

#primary embedding test of merged dataset suggested study as batch effect
#correct also for 2 sequencing batches of myelin mutant to avoid incomplete integration

#seperate batch 1 and batch2 from myelin mutant data before processing with integration pipeline
mic@meta.data$batch<-mic@meta.data$Study

mic@meta.data$batch[which(mic@meta.data$Study == "Myelin mutant")]<-mic@meta.data$SampleID[which(mic@meta.data$Study == "Myelin mutant")]

unique(mic@meta.data$batch)

mic@meta.data$batch<-factor(mic@meta.data$batch,
                           levels = unique(mic@meta.data$batch),
                           labels = c(rep("batch1",2), rep("batch2",2),"Zhou2020"))
unique(mic@meta.data$batch)

#prepare to correct for sample batch effect
#SCTranform pipeline

DefaultAssay(mic)<-"RNA"
seu.list<-SplitObject(mic, split.by = "batch")

for (k in c(1:length(seu.list))) {
    seu.list[[k]] <- SCTransform(seu.list[[k]], verbose = TRUE)#}

seu.features <- SelectIntegrationFeatures(object.list = seu.list, nfeatures = 500)
seu.list <- PrepSCTIntegration(object.list = seu.list, anchor.features = seu.features)

seu.anchors <- FindIntegrationAnchors(object.list = seu.list, normalization.method = "SCT", 
                                        anchor.features = seu.features, k.filter = 100)

seu.integrated <- IntegrateData(anchorset = seu.anchors, normalization.method = "SCT")

seu.integrated

#dimensionality reduction for integrated microglia dataset
seu.integrated <- RunPCA(seu.integrated, verbose = FALSE, npcs = 50)

ElbowPlot(seu.integrated, ndims = 50)

#neighbouring embedding for integrated dataset
seu.integrated <- RunUMAP(seu.integrated, dims = 1:20)

#visualize marker gene expression

microglia_marker<-c("Hexb","Cx3cr1","Aif1","Tmem119",
                   "C1qc","Trem2","Apoe","Cd68")

DefaultAssay(seu.integrated)<-"RNA"
seu.integrated<-NormalizeData(seu.integrated)
all.genes<-rownames(seu.integrated)
seu.integrated<-ScaleData(seu.integrated, features = all.genes)

for(k in microglia_marker){
    print(FeaturePlot(seu.integrated, feature = k))
}

#visualize integration result
#batch effect is removed and microglia subpopulations from individual dataset analysis
#is visually detectble at the embedded integrated dataset

DimPlot(seu.integrated, group.by = "Study")
DimPlot(seu.integrated, group.by = "Genotype")
DimPlot(seu.integrated, group.by = "CellType")
DimPlot(seu.integrated, group.by = "Age")
DimPlot(seu.integrated, group.by = "SampleID", label = TRUE)+NoLegend()

#to better confirm microglia subpopulations
#perform clustering analysis on integrated dataset

DefaultAssay(mic)<-"integrated"
mic<- FindNeighbors(object = mic, dims = 1:20)

mic <- FindClusters(object = mic, resolution = 0.2)
table(Idents(mic))

DimPlot(object = mic, reduction = 'umap',label = TRUE, 
        pt.size = 0.1,label.size = 6, repel = T) + theme(aspect.ratio=1)

mic@meta.data$CellType_redefine<-mic@meta.data$seurat_clusters

DefaultAssay(mic)<-"RNA"
Idents(mic)<-"seurat_clusters"

marker_all<-FindAllMarkers(mic, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.25,
                          test.use = "MAST")

#use DotPLot to confirm subpopulations

d2 <- by(marker_all, marker_all["cluster"], head, n=5)
d2<-as.matrix(d2)

d_sum2<-vector()
for (i in 1:length(d2)){
    d_temp<-as.character(d2[[i]]$gene)
    d_sum2<-append(d_sum2,d_temp)
}

d_sum2<-intersect(d_sum2,rownames(mic))

DotPlot(mic, features = d_sum2,
              dot.scale = 3.5
              #scale.by = "size"
             ) + coord_flip()+
theme(#strip.background = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 16,hjust = 0.5, face = "bold", angle = 90),
    axis.text.y = element_text(size = 8),
     legend.position = "right",
     #legend.spacing = unit(0, "mm"),
     legend.direction = "vertical",
        legend.text = element_text(size=5),
        legend.key.width = unit(4, "mm"),
        legend.key.height = unit(4, "mm"),
        legend.box.spacing = unit(1, "mm"),
        legend.margin = margin(2),
        legend.title = element_text(size = 7,angle = 90))

#check microglia subpopulation proportions in each genotype
Idents(mic)<-"seurat_clusters"

sep_pro<-as.data.frame(prop.table(table(mic@meta.data$seurat_clusters,
                                             mic@meta.data$Study_Genotype), margin = 2))

colnames(sep_pro)<-c("Cluster","Genotype","Proportion")


for(i in 0:5){
    test<-subset(sep_pro,subset = Cluster == i)
    test$Proportion<-round(test$Proportion,2)
    p<-ggplot(data=test, aes(x=Genotype, y=Proportion)) +
              geom_bar(stat="identity", fill=color[i+1])+
              geom_text(aes(label=Proportion), vjust=1.6, color="white", size=7)+
              theme_minimal()+
            ggtitle(paste0("Cluster ",i))+
    theme(plot.title = element_text(size=30),
         axis.text.x = element_text(size = 15,hjust = 0.5, angle = 90),
         axis.title.x = element_blank())
    print(p)
}

#rename clusters after marker gene confirmation
mic@meta.data$CellType_redefine<-factor(mic@meta.data$CellType_redefine,
                                       levels = c(0,1,2,3,4,5),
                                       labels = c("microglia1","DAM","microglia2",
                                                 "Myelin processing \n microglia",
                                                 "microglia1","unknown"))

######################visualization
#cell spatial location from each genotype
color<-c("#5445b1", "#749dae", "#f3c483", "#5c1a33", "#cd3341","#f7dc6a", "#4DA896", "#E07882")

sample_names<-unique(mic@meta.data$Study_Genotype)

sampleID_list<-vector(mode = "list", length = 6)

Idents(mic)<-"Study_Genotype"

for(i in 1:6){
    sampleID_list[[i]]<-WhichCells(mic, idents = sample_names[i])
}

for (k in 1:6){
       print(DimPlot(mic, label=F, group.by="Study_Genotype", 
        cells.highlight= sampleID_list[[k]], pt.size = 0.01)+ 
              scale_color_manual(labels = c("Other genotype", as.character(sample_names[k])), 
                                 values = c("grey", color[k])) +
              labs(color = "legend title"))
    }

#DimPlot of microglia subpopulations, color matched with manuscript
(plt<-DimPlot(mic, group.by = "CellType_redefine", label = T,
       cols = c("tan1","lightseagreen","darkgoldenrod2","lightskyblue", "grey48"))+NoLegend())



############################Peudotime trejactory calculation for DAM subpopulation
#subset DAM
seu.integrated<-subset(seu.integrated, subset = CellTYpe_redefine == "DAM")

DefaultAssay(seu.integrated)<-"integrated"
seu.sce<-as.SingleCellExperiment(seu.integrated)
seu.sce

#run trajectory calculation based on PCA embedding
#use prior DAM stage knowledge from Kerenschual et al 2017
#set microglia from Trem2KO_5xFAD as starting point (DAM1)
sds <- slingshot(seu.sce, reducedDim = "PCA", clusterLabels = seu.sce$Genotype, 
                 start.clus = "Trem2KO_5xFAD", 
                 allow.breaks = FALSE)

#visualize 
color<-c("#5445b1", "#749dae", "#f3c483", "#5c1a33",
         "#cd3341","#f7dc6a", "#4DA896", "#E07882")

plt + scale_fill_manual( values = color[c(5,2,6)])+
        theme_light()


















