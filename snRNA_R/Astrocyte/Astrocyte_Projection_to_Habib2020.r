#2021, Ting Sun
#myelin mutant snRNA-seq
#Astrocyte subset
#only CnpKO, MbpKO and WT
#differ from April 16 script:
    #object from myelin mutant was re-calculated after subset the genotypes
#query analysis with DAA (Habib 2020)
#Seurat workflow:https://satijalab.org/seurat/articles/integration_mapping.html

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

ast<-readRDS("./2021March_MyelinMutants_CnpKOandMbpKO_Astrocyte.rds")

daa<-readRDS( "./GSE143758_Admouse_Hippocampus_7m_AllNuclei_SCTint.rds")

#use Habib2020 dataset ask reference and query for cell subpopulation in myelin mutant
int.anchors<-FindTransferAnchors(reference = daa, query = ast, dims = 1:20)

predictions<-TransferData(anchorset = int.anchors, refdata = daa$orig.cluster, dims = 1:20)

ast.q<-AddMetaData(ast, metadata = predictions)

#check cell projection
#
ast.q

head(ast.q@meta.data)

#proportiuon tables
prop.table(table(ast.q@meta.data$Genotype,ast.q@meta.data$predicted.id), margin = 1)
prop.table(table(ast.q@meta.data$orig.cluster,ast.q@meta.data$predicted.id), margin = 1)
prop.table(table(ast.q@meta.data$orig.cluster,ast.q@meta.data$predicted.id), margin = 2)

#original paper DAA proportion
prop.table(table(daa@meta.data$Genotype,daa@meta.data$orig.cluster), margin = 1)

#spatial location
DimPlot(ast.q, group.by = "predicted.id")
DimPlot(ast.q, group.by = "orig.cluster")

#genotype UMAP location

color2<-c("orange","#749dae","#E07882")

sample_names<-unique(ast.q@meta.data$Genotype)
print(sample_names)

sampleID_list<-vector(mode = "list", length = length(sample_names))

Idents(ast.q)<-"Genotype"

for(i in 1:length(sample_names)){
    sampleID_list[[i]]<-WhichCells(ast.q, idents = sample_names[i])
}

for (k in 1:length(sample_names)){
        print(DimPlot(ast.q, label=F, group.by="Genotype", 
        cells.highlight= sampleID_list[[k]], pt.size = 0.01)+ 
              scale_color_manual(labels = c("Other Genotype", as.character(sample_names[k])), 
                                 values = c("grey", color2[k])) +
              labs(color = "legend title"))
    }

#precited cell proportion in myelin mutant
#Habib_4 is DAA

color<-c("grey48","peru","olivedrab3","turquoise3","steelblue3","deeppink1",
        "tan1","lightseagreen","darkgoldenrod2","lightskyblue", "tomato", "grey")

#export individual cell type proportion barplot
sep_pro<-as.data.frame(prop.table(table(ast.q@meta.data$predicted.id,
                                        ast.q@meta.data$Genotype), margin = 2))

colnames(sep_pro)<-c("CellType","Genotype","Proportion")

for(i in 1:6){
    test<-subset(sep_pro,subset = CellType == unique(sep_pro$CellType)[i])
    test$Proportion<-round(test$Proportion,3)
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



#export all required plots for manuscript
(plt1<-DimPlot(ast.q, group.by = "Genotype",
              cols = c("#749dae","#E07882","orange"))+
        labs(title = "Myelin Mutant AST Genotype UMAP"))

(predicted_daa<-DimPlot(ast.q, group.by = "predicted.id",
                        cols = c("gold2","tomato",
                                 "lightskyblue","lightseagreen", "coral4",
                                "rosybrown"))+
                labs(title = "DAA cluster prediction on MyelinMutant UMAP"))

mm_pro<-prop.table(table(ast.q@meta.data$Genotype,ast.q@meta.data$predicted.id), margin = 1)
ad_pro<-prop.table(table(daa@meta.data$Genotype,daa@meta.data$orig.cluster), margin = 1)

mm_pro
ad_pro








