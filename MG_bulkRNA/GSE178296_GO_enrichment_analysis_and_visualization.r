#Script for GSE178296
#MACS sorted microglia bulk RNA-seq resutls
#DEG gene ontology enrichment analysis
#and visualizations

#the script contains:
#enrichment analysis of pair comparisions using
    #ORA approach done by gprofiler2, GSEA was done by Webgestalt web interface
#aggregation and kmeans clustering of all DEGs from pair comparisons
    #gene list and clustering (based on scaled TPM)
    #enrichment of genes in each cluster using gprofiler2
    #complex heatmap visualization of enriched GO from all clusters

#run under r_backup_env


#attach packages

#for ORA enrichment
library(gprofiler2)

#for visualizations
library(dplyr)
library(ggplot2)
library(stringr)
library(viridis)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(cowplot)
library(simplifyEnrichment)

#reference database
library("org.Mm.eg.db")


#read in DEG result csv file lists
#previously applied cutoff of adjP<0.01

dir<-list.files(path = ".../indir/DESeq2_results/adjP<-0.01/",
               pattern = ".csv")

#generate matching name list
file_name<-c("5xFAD_vs_WT","Cnp_5xFAD_vs_5xFAD",
            "Cnp_5xFAD_vs_Cnp_ko", "Cnp_5xFAD_vs_WT",
            "Cnp_ko_vs_5xFAD", "Cnp_ko_vs_WT")

#create function and read in all files
filelist.read<-function(x){
    return(read.csv(file = paste0(".../indir/DESeq2_results/adjP<-0.01/",
                                  x), 
                    stringsAsFactors = FALSE))
}

file<-lapply(X = dir, FUN = ting.read)
names(file)<-file_name

#run for pathways (KEGG, WP and REAC)
for (i in 1:6){
    #set evcodes to True
    gostres_loop <- gost(query = file[[i]]$X, 
                organism = "mmusculus", ordered_query = FALSE, 
                multi_query = FALSE, significant = FALSE, exclude_iea = TRUE, 
                measure_underrepresentation = FALSE, evcodes = TRUE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = c("REAC","KEGG", "WP"), as_short_link = FALSE)

    
    #export result table
    write.csv(gostres_loop$result[,c(1:13,16)],
              file = paste0(".../outdir/out/",
                           file_name[i],"_GprofilerPathways_enrichment.csv"))
    
    #create GEM table and save
    gem <- gostres_loop$result[,c("term_id", "term_name", "p_value", "intersection")]
    colnames(gem) <- c("GO.ID", "Description", "p.Val", "Genes")
    gem$FDR <- gem$p.Val
    gem$Phenotype = paste0("+",i)
    gem <- gem[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Genes")]
    head(gem)
    
    write.table(gem,
             file = paste0(".../outdir/out/",
                          file_name[i],"_GprofilerPathways_GEM.txt"),
             sep = "\t", quote = F, row.names = F)
}

#run for GO:BP
for (i in 1:6){
    #set evcodes to True
    gostres_loop <- gost(query = file[[i]]$X, 
                organism = "mmusculus", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = TRUE, 
                measure_underrepresentation = FALSE, evcodes = TRUE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = c("GO:BP"), as_short_link = FALSE)
    
    #export result table
    write.csv(gostres_loop$result[,c(1:13,16)],
              file = paste0(".../outdir/out/",
                           file_name[i],"_GprofilerGOBP_enrichment.csv"))
    
    #create GEM table and save
    gem <- gostres_loop$result[,c("term_id", "term_name", "p_value", "intersection")]
    colnames(gem) <- c("GO.ID", "Description", "p.Val", "Genes")
    gem$FDR <- gem$p.Val
    gem$Phenotype = paste0("+",i)
    gem <- gem[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Genes")]
    head(gem)
    
    write.table(gem,
             file = paste0(".../outdir/out/",
                          file_name[i],"_GprofilerGOBP_GEM.txt"),
             sep = "\t", quote = F, row.names = F)
}

#pair comparison enrichment results were visualized using circular barplot refering to R plot gallery examples
#https://www.r-graph-gallery.com/web-circular-barplot-with-R-and-ggplot2.html


#for all DEGs, aggregate gene names with corresponding TPM values
#first read in all ranked DESeq2 results, cutoff applied: adjP<0.01
stat.list<-list.files(path = ".../indir/DESeq2_results/ranked_result")

setwd(".../indir/DESeq2_results/ranked_result")
stat<-lapply(stat.list,read.csv)
#note: all DESeq2 results have customized column names with comparison details (i.e. 5xFADvsWT_...)

#combine the list of DEG results into one
combined.stat<-stat.list %>% reduce(left_join, by = "ensembl_ID")



#read in TPM table for all genes across all sample replicates
tpm<-read.csv(".../indir/tpm.csv", stringsAsFactors = FALSE)
#read in corresponding meta data of sample replicates
#including original sequencing file name, replicate number, genotype and sample ID
meta.data<-read.csv(".../indir/meta_data/bulk_metadata.csv")

#filter tpm value based on significant genes targets (organized in combined.stat dataframe)
tpm<-subset(tpm, subset =  ENS_ID %in% combined.stat$ENS_ID)

#scale tpm matrix for further kmeans clustering analysis
x<-tpm #operate on a back up matrix and maintain the original one

x_norm<-as.data.frame(t(scale(t(x))))

#test suitable cluster number using sum of squeres
wss <- (nrow(x_norm)-1)*sum(apply(x_norm,2,var))
for (i in 1:40) wss[i] <- sum(kmeans(x_norm,centers=i)$withinss)
plot(1:40, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")

#perform kmeans clustering with k=10
#set reproducibility seed
set.seed(1001)
(cl <- kmeans(x_norm, 10))

#check cluster centers (with sample ID as rownames)
#and bind cluster center information with meta data
center<-t(cl$centers)
center.2<-cbind(meta.data,center)
#extract information of genes with matching cluster, each END_ID is assigned to one cluster
cluster<-data.frame(cl$cluster)
colnames(cluster)<-"adjp001_zcrore_k10"

#confirm filtered tpm matrix has identical rownames with cluster matrix (rownames based on ENS_ID)
all(rownames(tpm)==rownames(cluster))
#if TRUE
result<-cbind(tpm,cluster) #also works with merge based on rownames

#example script: use subeset function to extract one cluster
cluster1<-subset(result,subset = result$adjp001_zcrore_k10 == 1)

#export all cluster tables into a list
cluster_list<-vector(mode = "list", length = 10) #length = 10 for 10 clusters
for (i in 1:10){
    cluster_list[[i]]<-subset(result,subset = result$adjp001_zcrore_k10 == i)
}


#run each genes from each cluster for GO:BP
#save results
for (i in 1:10){
    #set evcodes to True
    gostres_loop <- gost(query = cluster_list[[i]]$Gene_symbol, 
                organism = "mmusculus", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = TRUE, 
                measure_underrepresentation = FALSE, evcodes = TRUE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = c("GO:BP"), as_short_link = FALSE)
    
    #export result table
    write.csv(gostres_loop$result[,c(1:13,16)],
              file = paste0(".../outdir/kmeans_gprofiler2/gprofiler_GOBP/outs/kmeans_",
                           file_name[i],"_GprofilerGOBP_enrichment.csv"))
    
    #create GEM table and save
    gem <- gostres_loop$result[,c("term_id", "term_name", "p_value", "intersection")]
    colnames(gem) <- c("GO.ID", "Description", "p.Val", "Genes")
    gem$FDR <- gem$p.Val
    gem$Phenotype = paste0("+",i)
    gem <- gem[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Genes")]
    head(gem)
    
    write.table(gem,
             file = paste0(".../outdir/kmeans_gprofiler2/gprofiler_GOBP/outs/kmeans_",
                          file_name[i],"GprofilerGOBP_GEM.txt"),
             sep = "\t", quote = F, row.names = F)
}



#aggregated enrichment analysis using GOBP enrichment result from each cluster
#re-read in GEM result from kmeans clusters
dir<-list.files(".../outdir/kmeans_gprofiler2/gprofiler_GOBP/outs/")

#prepare file names (caution, number system is 1, 10 ,2...9)
names<-paste0("cluster_",c(1,10,2:9))

#function for read in and organize result
read.go<-function(x){
    tmp<-read.csv(x, stringsAsFactors = FALSE)
    tmp<-subset(tmp, subset = p_value < 0.05)
    tmp<-tmp[,c("p_value","term_id","term_name")]
    tmp$sig_level<--log10(tmp$p_value)
    return(tmp)
}

setwd(".../outdir/kmeans_gprofiler2/gprofiler_GOBP/outs/")
data<-lapply(X = dir, FUN = read.go)
names(data)<-names

#final organization of colnames to suit for aggregate analysis
for(i in 1:10){
    colnames(data[[i]])[1]<-"p.adjust"
    colnames(data[[i]])[4]<-paste0(names[[i]],"_",colnames(data[[i]])[4])
    #data[[i]]<-data[[i]][,c(1,2)]
}

#reorder data
data<-data[c(1,3:10,2)]

#similarity and summary of GO terms from multiple GO results
p1<-simplifyGOFromMultipleLists(data, padj_cutoff = 0.05, ont = "BP",
                           db = "org.Mm.eg.db")
#color can be customized using heatmap_param argument






















