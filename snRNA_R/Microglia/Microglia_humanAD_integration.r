#2021, Ting Sun
#Myelin mutant scRAN-seq data
#Microglia subset

#merging with external datasets

#in house: Myelin mutant microglia profile

#mouse
#AD: Zhou 2020 7-month (5xFAD)

#human
#AD: Zhou 2020, Mathys 2019


###############after mathys 2019 dataset meta.data cleaning

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

#read in source data
mouse<-readRDS("./MyelinMutant_scRNAseq_Zhou_Microglia_raw_HUMANgene_TRANSLATED.rds")

#Zhou et.al 2020
zhou<-readRDS("./SYN21125841_Human_Zhou_MarcoColonna_NatMed_2020_AD_snRNA_microglia_scaled.rds")

#Mathys et.al 2019
mathys<-readRDS("./SYB18485175_Human_Mathys_Tsai_Nature_2019_AD_snRNA_filter_scaled.rds")

#meta.data organization and clean up
#remaining tags: 
#orig.ident, nCount_RNA, nFeature_RNA, percent.mt, CellType, SampleID, Genotype, Age, Tissue, Study, Sex


#dataset merging
obj<-list(mouse, zhou, mathys)

mic<-merge(x = obj[[1]], y = c(obj[2:3]), project = "Depp_Sun2021")

mic
head(mic@meta.data)

#check all meta.data entries

for(i in c(4,5,7,8,9,10)){
    print(unique(mic@meta.data[,i]))
}

table(mic@meta.data$Study, mic@meta.data$Species)

mic@meta.data$Condition[which(mic@meta.data$Condition =="AD")]<-paste0("AD_",
                                                                       mic@meta.data$Study[which(mic@meta.data$Condition=="AD")])

mic@meta.data$Condition[which(mic@meta.data$Condition =="Ctrl")]<-paste0("Ctrl_",
                                                                         mic@meta.data$Study[which(mic@meta.data$Condition=="Ctrl")])

#primary embedding analysis
#detect potential batch effect
#one step normalisation until PCA embedding
DefaultAssay(mic)<-"RNA"

mic<-NormalizeData(mic)

mic<-FindVariableFeatures(mic, selection.method = 'vst',
                            nfeatures = 2000)

all.genes <- rownames(mic)
mic <- ScaleData(mic, features = all.genes)

mic<-RunPCA(mic, features = VariableFeatures(object = mic))
ElbowPlot(mic, ndims = 50)

mic<-RunUMAP(mic, dims = 1:30)

DimPlot(mic, group.by = "Study")
DimPlot(mic, group.by = "Condition")
DimPlot(mic, group.by = "Species")
DimPlot(mic, group.by = "SampleID", label = T)+NoLegend()

#preare to correct study and sequnecing batch effects

#seperate batch 1 and batch2 from myelin mutant data before processing with integration pipeline
mic@meta.data$batch<-mic@meta.data$Study

mic@meta.data$batch[which(mic@meta.data$Study == "Myelin mutant")]<-mic@meta.data$SampleID[which(mic@meta.data$Study == "Myelin mutant")]

unique(mic@meta.data$batch)

mic@meta.data$batch<-factor(mic@meta.data$batch,
                           levels = unique(mic@meta.data$batch),
                           labels = c(rep("batch1",4), rep("batch2",4),"Zhou2020","Mathys2019"))
unique(mic@meta.data$batch)

mic@meta.data$batch<-paste0(mic@meta.data$batch,mic@meta.data$Species)

unique(mic@meta.data$batch)

table(mic@meta.data$batch)

#prepare to correct for sample batch effect
#SCTranform pipeline

DefaultAssay(mic)<-"RNA"
seu.list<-SplitObject(mic, split.by = "batch")

for (k in c(1:length(seu.list))) {
    seu.list[[k]] <- SCTransform(seu.list[[k]], verbose = TRUE)
}

seu.features <- SelectIntegrationFeatures(object.list = seu.list, nfeatures = 500)
seu.list <- PrepSCTIntegration(object.list = seu.list, anchor.features = seu.features)

seu.anchors <- FindIntegrationAnchors(object.list = seu.list, normalization.method = "SCT", 
                                        anchor.features = seu.features, 
                                      k.filter = 300)

seu.integrated <- IntegrateData(anchorset = seu.anchors, normalization.method = "SCT")

seu.integrated

seu.integrated <- RunPCA(seu.integrated, verbose = FALSE, npcs = 50)

ElbowPlot(seu.integrated, ndims = 50)

seu.integrated <- RunUMAP(seu.integrated, dims = 1:15)

#check for marker gene expression and removal of batch effect

microglia_marker<-c("Hexb","Cx3cr1","Aif1","Tmem119",
                   "C1qc","Trem2","Apoe","Cd68")

microglia_marker<-toupper(microglia_marker)

DefaultAssay(seu.integrated)<-"RNA"
seu.integrated<-NormalizeData(seu.integrated)
all.genes<-rownames(seu.integrated)
seu.integrated<-ScaleData(seu.integrated, features = all.genes)

for(k in microglia_marker){
    print(FeaturePlot(seu.integrated, feature = k))
}

DimPlot(seu.integrated, group.by = "Study")
DimPlot(seu.integrated, group.by = "Condition")
DimPlot(seu.integrated, group.by = "CellType")
DimPlot(seu.integrated, group.by = "Age")
DimPlot(seu.integrated, group.by = "SampleID", label = TRUE)+NoLegend()

#perform clustering

DefaultAssay(seu.integrated)<-"integrated"
seu.integrated<- FindNeighbors(object = seu.integrated, dims = 1:15)

seu.integrated <- FindClusters(object = seu.integrated, resolution = 0.5)
table(Idents(seu.integrated))

DimPlot(object = seu.integrated, reduction = 'umap',label = TRUE, 
        pt.size = 0.1,label.size = 6, repel = T) + theme(aspect.ratio=1)

seu.integrated@meta.data$seurat_clusters<-seu.integrated@meta.data$integrated_snn_res.0.5

#marker gene calculation
DefaultAssay(seu.integrated)<-"RNA"
Idents(seu.integrated)<-"seurat_clusters"

marker_all<-FindAllMarkers(seu.integrated, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.25,
                          test.use = "MAST")

#detect and plot cell subpopulation proportion changes of each genotype
#visualize top5 to 10 marker genes in DopPlot
d2 <- by(marker_all, marker_all["cluster"], head, n=5)
d2<-as.matrix(d2)

d_sum2<-vector()
for (i in 1:length(d2)){
    d_temp<-as.character(d2[[i]]$gene)
    d_sum2<-append(d_sum2,d_temp)
}

d_sum2<-intersect(d_sum2,rownames(seu.integrated))

DotPlot(seu.integrated, features = d_sum2,
              dot.scale = 3.5
              #scale.by = "size"
             ) + coord_flip()+
theme(#strip.background = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 10,hjust = 0.5, face = "bold"),
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

#very rough proportion analysis for DAM, prunning and myelin processing microglia
seu.integrated@meta.data$temp<-seu.integrated@meta.data$seurat_clusters

seu.integrated@meta.data$temp<-factor(seu.integrated@meta.data$temp,
                                     levels = c("0","2","3","6","7","8","9","10","1","4","5"),
                                     labels = c(rep("other microglia",8),
                                                "DAM", "Prunning","Myelin processing"))

DimPlot(seu.integrated, group.by = "temp")

unique(seu.integrated@meta.data$Condition)

#calculate subpopulation shift based on cluster
Idents(seu.integrated)<-"temp"

sep_pro<-as.data.frame(prop.table(table(seu.integrated@meta.data$temp,
                                             seu.integrated@meta.data$Condition), margin = 2))
head(sep_pro)

colnames(sep_pro)<-c("CellType", "Condition", "Proportion")

unique(sep_pro$CellType)
unique(sep_pro$Condition)

color<-c("tomato","springgreen4","turquoise3","magenta")

sep_pro$Condition<-factor(sep_pro$Condition,
                         levels = c("WT","Ctrl_Mathys2019","Ctrl_Zhou2020","AD_Mathys2019","AD_Zhou2020",
                                    "AD_TREM2_R62H","5xFAD","Trem2KO_5xFAD","Trem2KO",
                                   "CnpKO","Plp1KO","Foxg1MbpKO"))

sep_pro2<-subset(sep_pro, subset = Condition %in% c("WT","Ctrl_Mathys2019","Ctrl_Zhou2020","AD_Mathys2019","AD_Zhou2020",
                                                    "AD_TREM2_R62H","5xFAD",
                                   "CnpKO","Plp1KO","Foxg1MbpKO"))

for(i in 1:4){
    test<-subset(sep_pro2,subset = CellType == unique(sep_pro$CellType)[i])
    test$Proportion<-round(test$Proportion,2)
    p<-ggplot(data=test, aes(x=Condition, y=Proportion)) +
              geom_bar(stat="identity", fill=color[i])+
              geom_text(aes(label=Proportion), vjust=1.6, color="white", size=7)+
              theme_minimal()+
            ggtitle(unique(sep_pro$CellType)[i])+
    theme(plot.title = element_text(size=30),
         axis.text.x = element_text(size = 15,hjust = 0.5, angle = 90),
         axis.title.x = element_blank())
    print(p)
}

########################color code matched with manuscript

(plt<-DimPlot(seu.integrated, group.by = "temp", label = T, label.size = 6, repel = T,
       cols = c("tan1","lightseagreen","rosybrown3","lightskyblue"))+
NoLegend()+
ggtitle("Microglia human mouse merged object, new color"))

#highlight genotype with manuscript matched color
color<-c("#5445b1", "#749dae", "#cd3341",
         "#5c1a33","#f3c483",  "#4DA896", "#E07882","#f7dc6a")

sample_names<-c("WT","CnpKO","5xFAD",
                "AD_Zhou2020","Ctrl_Zhou2020",
                "AD_Mathys2019","Ctrl_Mathys2019")

sampleID_list<-vector(mode = "list", length = 7)

Idents(seu.integrated)<-"Condition"

for(i in 1:7){
    sampleID_list[[i]]<-WhichCells(seu.integrated, idents = sample_names[i])
}

for (k in 1:7){
        print(DimPlot(seu.integrated, label=F, group.by="Condition", 
        cells.highlight= sampleID_list[[k]])+ 
              scale_color_manual(labels = c("Other Conditions", as.character(sample_names[k])), 
                                 values = c("grey", color[k])) +
              labs(color = "legend title"))
    }








