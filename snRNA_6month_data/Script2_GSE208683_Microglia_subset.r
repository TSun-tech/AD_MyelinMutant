#GSE208683, analysis script 2
#data collected from snRNA-seq of mouse brain hemispheres
#in total 4 genotypes: WT, CnpKO, 5xFAD, CnpKO 5xFAD

#this script contains analysis of microglia subset from the data
#in which the microglia subpopulations were determined,
#and calculated for marker genes, proportions across genotypes

#run under sc_update, R 4.1.2, Seurat 4.1.1

#attach analysis packages

#for scRNA-seq objects
library(Seurat)
library(dplyr)
library(Matrix)
library(abind)
library(sctransform)
library(SingleCellExperiment)

#for visualization
library(cowplot)
library(ggplot2)

#for DGE analysis
library(MAST)
library(patchwork)
#library(pheatmap)

#read in processed object from all cells
seu<-readRDS(".../indir/scRNAseq_clustered.rds")

#subset microglia and macrophages for zoom in analysis
mic<-subset(seu, subset = CellType %in% c("Microglia", "Macrophages"))

#apply SCTransform normalizion before dimentionality reduction
DefaultAssay(mic)<-"RNA"

mic<-SCTransform(mic)

#linear dimention reduction using PCA
mic<-RunPCA(mic)
ElbowPlot(mic, ndims = 50)

#based on elbowplot, use first 20PC for neighbouring embedding
#visualize UMAP
mic<-RunUMAP(mic, dims = 1:20)

DimPlot(mic, group.by = "Genotype")
DimPlot(mic, group.by = "CellType")
DimPlot(mic, group.by = "seurat_clusters", label = T)+NoLegend()

#perform unbiased clustering and calculate cluster marker genes
#detect potential noise clusters
mic <- FindNeighbors(mic, dims = 1:20)

mic<-FindClusters(mic, resolution = 0.5)

table(Idents(mic))

DimPlot(object = mic, reduction = 'umap',label = TRUE, 
        pt.size = 0.1,label.size = 6, repel = T) + theme(aspect.ratio=1)+NoLegend()

Idents(mic)<-"seurat_clusters"

markers_sct<-FindAllMarkers(mic, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0)

#in total 8 clusters were detected under resolution 0.5
#cluster 6 to 7 are identified as doublets or contaminated cells
#clean dataset and perform reanalysis
mic_clean<-subset(mic, subset = seurat_clusters %in% c(0:5))

mic_clean<-SCTransform(mic_clean)

#proceed woth embedding
mic_clean<-RunPCA(mic_clean)
ElbowPlot(mic_clean, ndims = 50)

mic_clean<-RunUMAP(mic_clean, dims = 1:40) #higher amount of PC that is recommended by SCTransform pipeline

DimPlot(mic_clean, group.by = "Genotype")
DimPlot(mic_clean, group.by = "CellType")
DimPlot(mic_clean, group.by = "seurat_clusters", label = T)+NoLegend()

#annotate seurat clusters based on calculated marker genes
mic@meta.data$Cell_subtype<-factor(mic@meta.data$Cell_subtype,
                                  levels = c(0:5),
                                  labels = c("AD_DAM",
                                             "Homeostatic", "Homeostatic",
                                            "MyTE", 
                                             "Myelin_DAM", 
                                             "Macrophages"))

#calculate subpopulation quantity (absolute quantity) across each genotype
proportion<-prop.table(table(mic@meta.data$Genotype, mic@meta.data$Cell_subtype), margin = 1)

abs<-table(mic@meta.data$SampleID, mic@meta.data$Cell_subtype)

#DE analysis between subpopulation and identities
#example script using Myelin DAM vs AD DAM
#test threshold only using min.pct 0.1
Idents(mic)<-"Cell_subtype"
dam_deg<-FindMarkers(mic, ident.1 = "Myelin_DAM", ident.2 = "AD_DAM", min.pct = 0.1, logfc.threshold = 0)

#example script for visualizing DEG result using volcano plot
# The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)

# add a column of NAs
dam_deg$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
dam_deg$diffexpressed[dam_deg$avg_log2FC > 0 & dam_deg$p_val_adj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
dam_deg$diffexpressed[dam_deg$avg_log2FC < 0 & dam_deg$p_val_adj < 0.05] <- "DOWN"

# Re-plot but this time color the points with "diffexpressed"
p <- ggplot(data=dam_deg, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed)) + geom_point() + theme_minimal()

# Add lines that indicates threshold
(p2 <- p + geom_vline(xintercept=0, col="grey48", linetype = "longdash") +
        geom_hline(yintercept=-log10(0.05), col="red"))

# Now write down the name of genes besidam_deg the points...
# Create a new column "dam_deglabel" to dam_deg, that will contain the name of genes differentially expressed (NA in case they are not)
dam_deg$dam_deglabel <- NA
dam_deg$dam_deglabel[dam_deg$diffexpressed != "NO"] <- dam_deg$X[dam_deg$diffexpressed != "NO"]

ggplot(data=dam_deg, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=dam_deglabel)) + 
    geom_point() + 
    theme_minimal() +
    geom_text_repel(max.overlaps = 20) +
    scale_color_manual(values=c("lightblue3", "grey", "lightcoral"))+ 
geom_vline(xintercept=0, col="grey48", linetype = "longdash") +
        geom_hline(yintercept=-log10(0.05), col="grey48", linetype = "longdash")

#selective labeling of genes
dam_deg$dam_deglabel <- NA
position<-c(dam_deg$X %in% c(#upregulate AD DAM
    "Plxdc2","Csmd3","Nav3","Ank","Myo1e","Hdac9","Cacna1a","Ctnna3","Apbb1ip","Gnas",
                            "Dgkd","Sdk1","Serpine2",
                            #upregulate myelinDAM
    "Apobec1","Gpnmb","Ms4a7","Lyz2","Apoe", "Abca1", "Cd84", "Pid1", "Fyb", "Ms4a6d", "Fyb",
    "Mdfic","Ctsh","Ccnd3")
)
dam_deg$dam_deglabel[position] <- dam_deg$X[position]

#plot
(p_final<-ggplot(data=dam_deg, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=dam_deglabel)) + 
    geom_point() + 
    theme_minimal() +
    xlim(c(-1.2,1.2))+
    geom_text_repel(max.overlaps = 20, size=6.5, face = "bold") +
    scale_color_manual(values=c("plum4", "grey", "deepskyblue3"))+ 
    geom_vline(xintercept=0, col="grey48", linetype = "longdash") +
        geom_hline(yintercept=-log10(0.05), col="grey48", linetype = "longdash")
 )












