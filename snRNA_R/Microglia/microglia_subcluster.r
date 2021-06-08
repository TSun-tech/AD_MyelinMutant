#2021April06, Ting Sun
#myelin mutants all samples
#Microglia subset and DGE analysis
#CnpKO, Mbp cKO and corresponding WT

#to show enrichment of myelin processing microglia in both CnpKO and Mbp cKO
#in Mbp cKO to greater extend
#run under sc_env

#attach analysis packages
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

#read in all cells and prepare subset
#seu<-readRDS(./2020Oct_MyelinMutants_all_cleaned_run2.rds")

#subset microglia
mic<-subset(seu, subset = CellType %in% c("Microglia") & Genotype %in% c("CnpKO", "WT"))

#proceed with microglia concentrated analysis
#using 2000 varialble features

DefaultAssay(mic)<-"RNA"

mic<-NormalizeData(mic)

mic<-FindVariableFeatures(mic, selection.method = 'vst',
                            nfeatures = 2000)

all.genes<-rownames(mic)

mic<-ScaleData(mic, features = all.genes)

#after normalization and scaling, validate KO genes in individual replicates

for (gene in c("Mbp","Cnp","Plp1")){
    print(VlnPlot(mic, features = gene, group.by = "SampleID"))
    print(RidgePlot(mic, features = gene, group.by = "SampleID"))
}

#proceed woth embedding
mic<-RunPCA(mic, features = VariableFeatures(object = mic))
ElbowPlot(mic, ndims = 50)

#select first 20 PCs for microglia neighbouring embedding
mic<-RunUMAP(mic, dims = 1:15, #set.op.mix.ratio = 0.1, 
              min.dist = 0.3, spread = 0.3, 
              verbose = FALSE)

#check embedding with kockout genes
DimPlot(mic, group.by = "Genotype")
DimPlot(mic, group.by = "CellType")

#evaluate microglia cell annotation by highlighting the marker genes
microglia_marker<-c("Hexb","Cx3cr1","Aif1","Tmem119",
                   "C1qc","Trem2","Apoe","Cd68","Il23a","Il12a","Il12b")

for(k in microglia_marker){
    print(FeaturePlot(mic, feature = k))
}

#confirm that the spatial seperation was not caused by 3 kockout genes
#extract raw count matrix
exp<-as.matrix(mic@assays$RNA@counts)

p1<-match("Cnp",rownames(exp))

exp_noko<-exp[-p1,]

mic_new<-CreateSeuratObject(exp_noko, meta.data= mic@meta.data, 
                       assay = "RNA")

DefaultAssay(mic_new)<-"RNA"

mic_new<-NormalizeData(mic_new)

mic_new<-FindVariableFeatures(mic_new, selection.method = 'vst',
                            nfeatures = 2000)

all.genes<-rownames(mic_new)

mic_new<-ScaleData(mic_new, features = all.genes)

mic_new<-RunPCA(mic_new, features = VariableFeatures(object = mic_new))
ElbowPlot(mic_new, ndims = 50)

mic_new<-RunUMAP(mic_new, dims = 1:15, #set.op.mix.ratio = 0.1, 
              min.dist = 0.3, spread = 0.3, 
              verbose = FALSE)

#embedding without 3 ko genes
DimPlot(mic_new, group.by = "Genotype")
DimPlot(mic_new, group.by = "CellType")
###########result suggested cell neighbouring embedding is not related with KO genes

##################################
#after checking embedding without KO genes
#proceed with clustering analysis and cell subpopulation annotation

#find neibours
#set.seed(20201217)

mic <- FindNeighbors(object = mic, dims = 1:15)

#test runs suggested clustering at resolution 1.5

mic <- FindClusters(object = mic, resolution = 1.5)
table(Idents(mic))

DimPlot(object = mic, reduction = 'umap',label = TRUE, 
        pt.size = 0.1,label.size = 6, repel = T) + theme(aspect.ratio=1)

#calculate subcluster marker genes
Idents(mic)<-"RNA_snn_res.1.5"
DefaultAssay(mic)<-"RNA"

markers<-FindAllMarkers(mic,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, 
                        test.use = "MAST")

#annotate microglia subclusters based on calculated marker gene list
mic@meta.data$CellType<-mic@meta.data$seurat_clusters

mic@meta.data$CellType<-factor(mic@meta.data$CellType,
                              levels = c(0,1,2,3,4,5,6,7,8),
                              labels = c("Homeostatic", "Homeostatic","Homeostatic",
                                        "Myelin processing","Homeostatic","Homeostatic",
                                         "DAM","Brain Border \n Marcrophages",
                                        "unknown"))

DimPlot(mic, group.by = "CellType", label = T)+NoLegend()



#calculate relative proportion of microglia subpopulations in each genotype
#plot individual barplot for each population

(per<-prop.table(table(mic@meta.data$Genotype, mic@meta.data$CellType), margin = 1))

Idents(mic)<-"CellType"

#export individual cell type proportion barplot
sep_pro<-as.data.frame(prop.table(table(mic@meta.data$CellType,mic@meta.data$Genotype), margin = 2))

colnames(sep_pro)<-c("CellType","Genotype","Proportion")

color<-c("tan1","lightskyblue","lightcoral","lightseagreen","grey48")

for(i in 1:5){
    test<-subset(sep_pro,subset = CellType == unique(sep_pro$CellType)[i])
    test$Proportion<-round(test$Proportion,2)
    p<-ggplot(data=test, aes(x=Genotype, y=Proportion)) +
              geom_bar(stat="identity", fill=color[i])+
              geom_text(aes(label=Proportion), vjust=1.6, color="white", size=7)+
              theme_minimal()+
            ggtitle(unique(sep_pro$CellType)[i])+
    theme(plot.title = element_text(size=30),
         axis.text.x = element_text(size = 20,hjust = 0.5),
         axis.title.x = element_blank())
    print(p)
}


################backup code 1
################color scheme matched with publication
color_oder <- setNames(
    c("tan1","lightskyblue","lightcoral","lightseagreen","grey48"),
    c("Homeostatic", "Myelin processing","DAM",
      "Brain Border \n Marcrophages","unknown"))

DimPlot(mic, group.by = "CellType",
        cols = c("tan1","lightskyblue","lightseagreen","lightcoral","grey48"),
       label = T, repel = T)+NoLegend()

################backup code 2
################DimPlot highlighting cells from each genotye/replicate

sample_names<-unique(mic@meta.data$Genotype)

sample_names

sample_names<-factor(sample_names, levels = unique(sample_names))

sampleID_list<-vector(mode = "list", length = 2)

Idents(mic)<-"Genotype"

for(i in 1:2){
    sampleID_list[[i]]<-WhichCells(mic, idents = sample_names[i])
}

color<-c("#5445b1", "#749dae", "#E07882")

for (k in 1:2){
    print(DimPlot(mic, label=F, group.by="Genotype", 
        cells.highlight= sampleID_list[[k]], pt.size = 0.01)+ 
              scale_color_manual(labels = c("Other genotype", as.character(sample_names[k])), 
                                 values = c("grey", color[k])) +
              labs(color = "legend title"))
    }








