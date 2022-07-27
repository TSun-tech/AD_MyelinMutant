#Depp and Sun et al 2022, analysis script 3
#reanalysis of human snRNAseq data microglia subsets

#data collected from Mathys et al 2019 (SYB18485175),
#and Zhou et al 2020(SYN21670836)

#script contains peremetes and filter that used for microglia analysis from both datasets

#run under sc_update, R 4.1.2, Seurat 4.1.1

#attach packages
library(Seurat)
library(dplyr)
library(ggplot2)

#read in microglia subsets from both human snRNAseq datasets
zhou<-readRDS(".../indir/rds/human/SYN21125841_Human_Zhou_MarcoColonna_NatMed_2020_AD_snRNA_microglia_IntReady.rds")
mathys<-readRDS(".../indir/rds/SYB18485175_Human_Mathys_Tsai_Nature_2019_AD_snRNA_onlyMG_IntReady.rds")

#proceed with Mathys et al 2019 dataset first

#observe sample qualities
FeatureScatter(mathys, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "SampleID")

#proceed dataset with standard analysis until PCA
#define suitable nPCs and run clustering analysis
#calculate cluster markers and determine potential noise cell populations

DefaultAssay(mathys)<-"RNA"

mathys<-NormalizeData(mathys)

mathys<-FindVariableFeatures(mathys, selection.method = 'vst',
                            nfeatures = 2000)

all.genes <- rownames(mathys)
mathys <- ScaleData(mathys, features = all.genes)

mathys<-RunPCA(mathys, features = VariableFeatures(object = mathys))
ElbowPlot(mathys, ndims = 50)

#run neighbouring ebedding
mathys<-RunUMAP(mathys, dims = 1:40, verbose = FALSE)

#continue obstain subclusters and marker genes
mathys<-FindNeighbors(mathys, dims = 1:40)

mathys<-FindClusters(mathys, resolution = 0.5)

DimPlot(object = mathys, reduction = 'umap',label = TRUE, 
        pt.size = 0.1,label.size = 6, repel = T) + theme(aspect.ratio=1)+NoLegend()

Idents(mathys)<-"seurat_clusters"

markers<-FindAllMarkers(mathys, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0)

#remove identified doublet and noise clusters
#and annotated remaining cell subpopulations
mathys<-subset(mathys, subset = seurat_clusters %in% c(0,1,2,3,4))

mathys@meta.data$Cell_subtype<-mathys@meta.data$seurat_clusters

mathys@meta.data$Cell_subtype<-factor(mathys@meta.data$Cell_subtype,
                                      levels = c(1,0,2,4,3),
                                     labels = c("Hom1","Hom2","DAM",
                                               "Ribo_high","Leukocyte"))

#function for checking gene expressions on UMAP
featureplot_pub<-function(x){
    print(FeaturePlot(object = mathys,features = x,
        reduction = 'umap',label = F)+
          theme(axis.line=element_blank(), axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position = "none",
          title = element_text(size = 50, face = "bold"))+
scale_color_gradient(low = "lemonchiffon", high = "skyblue4")+ 
            scale_fill_continuous(limits=c(0, 6), breaks=seq(0,6,by=7)))
}

#check co-expressions of two genes
#example script using MS4A6A and SPP1
(combine_plot<-FeaturePlot(mathys, features = c("MS4A6A", "SPP1"), blend = TRUE, blend.threshold = 0, pt.size = 2,
           order = TRUE, combine = FALSE,
            #cols = c("blue", "yellow", "green")
           )#+NoAxes()
 )



####################################################
#second part od the script focus on reanalysis of Zhou et al 2020 snRNAseq data microglia subset

#check meta data labeling
head(zhou@meta.data.data)

#check data quality
zhou[["percent.mt"]]<- PercentageFeatureSet(zhou, pattern = "^MT-")

Idents(zhou)<-"Condition"
VlnPlot(zhou, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(zhou, feature1 = "nCount_RNA", feature2 = "percent.mt")

#apply one step additional filter
zhou <- subset(zhou, subset = nFeature_RNA > 200 & nCount_RNA < 10000 & percent.mt < 20) 

#select only sporadic AD cases and remove identifies outlier sample
zhou<-subset(zhou, subset = SampleID != "C1" & Condition != "AD_TREM2_R62H")

#proceed with data normalization, dimentionality reduction and neighbouring embedding
DefaultAssay(zhou)<-"RNA"

zhou<-NormalizeData(zhou)

zhou<-FindVariableFeatures(zhou, selection.method = 'vst',
                            nfeatures = 2000)

all.genes <- rownames(zhou)
zhou <- ScaleData(zhou, features = all.genes, vars.to.regress = c("nCount_RNA", "nFeature_RNA","percent.mt"))

zhou<-RunPCA(zhou, features = VariableFeatures(zhou))
ElbowPlot(zhou, ndims = 50)

zhou<-RunUMAP(zhou, dims = 1:30, verbose = FALSE)

#perform unbiased clustering and remove noise clusters

zhou<-FindNeighbors(zhou, dims = 1:30)

zhou<-FindClusters(zhou, resolution = 0.5)

DimPlot(object = zhou, reduction = 'umap',label = TRUE, 
        pt.size = 0.1,label.size = 6, repel = T) + theme(aspect.ratio=1)+NoLegend()

#remove cluster 5 (ribo high), and cluster 6 (mito high)
#futher reanalyze cleaned dataset

zhou<-subset(zhou, subset = seurat_clusters %in% c(5,6,7), invert = TRUE)

DefaultAssay(zhou)<-"RNA"

zhou<-NormalizeData(zhou)

zhou<-FindVariableFeatures(zhou, selection.method = 'vst',
                            nfeatures = 2000)

all.genes <- rownames(zhou)
zhou <- ScaleData(zhou, features = all.genes, vars.to.regress = c("nCount_RNA", "nFeature_RNA","percent.mt"))

zhou<-RunPCA(zhou, features = VariableFeatures(zhou))
ElbowPlot(zhou, ndims = 50)

zhou<-RunUMAP(zhou, dims = 1:30, verbose = F)


zhou@meta.data$Cell_subtype<-zhou@meta.data$seurat_clusters
zhou@meta.data$Cell_subtype<-factor(zhou@meta.data$Cell_subtype,
                                   levels = c(1,0,2,3,4,9,8),
                                   labels = c("Hom1","Hom2","DAM","Histone_high",
                                             "MyTE", "Monocyte",
                                              "Leukocyte"))

#visualize condition distributions on UMAP
zhou@meta.data$Condition<-factor(zhou@meta.data$Condition,
                                levels = c("Ctrl","AD"))

(plt1<-DimPlot(zhou, group.by = "Condition", cols = c("gray70", "plum4"),
                pt.size = 1.2, label.size = 5, label = F, repel = T)+
          theme(axis.line=element_blank(), axis.text.x=element_blank(),
              axis.text.y=element_blank(),axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),legend.position = "right",
              title = element_text(size = 20, face = "bold"))+
          ggtitle("Zhou Microglia Conditions")
)
#similar to Mathys dataset, export single and combined gene expressions on UMAP

#function for single gene expression in Zhou 2020 dataset
featureplot_pub<-function(x){
    print(FeaturePlot(object = zhou,features = x,
        reduction = 'umap',label = F)+
          theme(axis.line=element_blank(), axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position = "none",
          title = element_text(size = 50, face = "bold"))+
scale_color_gradient(low = "lemonchiffon", high = "skyblue4")+ 
            scale_fill_continuous(limits=c(0, 6), breaks=seq(0,6,by=7)))
}

    #combined gene expressions
(combine_plot<-FeaturePlot(zhou, features = c("MS4A7", "TREM2"), blend = TRUE, blend.threshold = 0, pt.size = 2,
           order = TRUE, combine = FALSE,
            #cols = c("blue", "yellow", "green")
           )#+NoAxes()
 )


















