#2020May06, Ting Sun
#AD microglia bulk RNAseq data
#step 2 for kmeans clustering and heatmap visualization
#after tested clustering using log2 transformed TPM value (see script)
#result still showed high variants on average expression level and clustering were influenced
#try with z-score normalised matrix

#reference: https://www.r-bloggers.com/k-mean-clustering-heatmap/
#   https://github.com/jokergoo/ComplexHeatmap/issues/136

#analysis under visual env

library(dplyr)
library(colorspace) 
library(dendextend)
library(circlize)
library(gplots)
library(ggplot2)
library(corrplot)
library(pheatmap)

#read in prepared TPM matrix (pre-filtered with adjP < 0.01) 
#and meta.data (for plot), match sample ID with genotypes
df<-read.csv("./2020_AD_microglia_TPM_adjp001_for_kmenas.csv")
meta.data<-read.csv("./bulk_metadata.csv")

head(df)
head(meta.data)
dim(df)

rownames(df)<-df$ENS_ID

x<-df[,6:21]

colnames(x)<-gsub("X5xFAD_","5xFAD_",colnames(x))
colnames(x)<-gsub("_TPM","",colnames(x))

head(x)

#normalise data with z-score
x_norm<-as.data.frame(t(scale(t(x))))
head(x_norm)
barplot(x_norm$WT_1)

#test suitable cluster number using sum of squeres
wss <- (nrow(x_norm)-1)*sum(apply(x_norm,2,var))
for (i in 1:40) wss[i] <- sum(kmeans(x_norm,centers=i)$withinss)
plot(1:40, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")

wss

hc <- hclust(dist(x), "ward")
plot(hc) # the plot can also help to decide the # of clusters
memb <- cutree(hc, k = 10)
head(memb,30)

#perform kmeans clusering with k = 10
#extract cluster center for each replicate

#set.seed(16593)
(cl <- kmeans(x_norm, 10))

str(cl)

cl$centers

center<-t(cl$centers)
head(center)

#further refine matrix for plotting purpose

colnames(center)<-paste0("cluster_",colnames(center))
head(center)

all(meta.data$sample==rownames(center))

center.2<-cbind(meta.data,center)

center.2

colnames(center.2)[6:15]

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
     cluster = col)
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  #data_sum <- rename(data_sum, c("mean" = varname))
 return(data_sum)
}

empty.list<-vector(mode = "list", length = 10)

empty.list

(df2<-data_summary(center.2, varname="cluster_1", 
                    groupnames=c("genotype")))

(df3<-data_summary(center.2, varname="cluster_2", 
                    groupnames=c("genotype")))

(df_combine<-rbind(df2,df3))

for (i in 1:10){
    empty.list[[i]]<-data_summary(center.2, varname= colnames(center.2)[6:15][i], 
                    groupnames=c("genotype"))
}

empty.list

fin.table<-rbind_list(empty.list)

fin.table

fin.table$mean<-as.numeric(fin.table$mean)
fin.table$sd<-as.numeric(fin.table$sd)
head(fin.table)

#define plot factor levels

fin.table$cluster<-factor(fin.table$cluster,
                         levels = c("cluster_1", "cluster_2", "cluster_3", "cluster_4", "cluster_5",
                                   "cluster_6","cluster_7","cluster_8","cluster_9","cluster_10"))
fin.table$genotype<-factor(fin.table$genotype,
                          levels = c("WT", "5xFAD", "Cnp_ko","Cnp_5xFAD"))


#plot cluster centers with SD calculated with replicates
p<- ggplot(fin.table, aes(x=genotype, y=mean, group=cluster, color=cluster)) + 
  geom_line() +
  geom_point()+
  ylab("Average center value")+
    theme(text = element_text(size=20))+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                 position=position_dodge(0.05))+
    scale_color_manual(values = c("grey50","coral1","yellowgreen","yellow3","maroon3",
                                "cadetblue","brown","lightpink","steelblue4","black"))
print(p)

#export line plot seperately

color_vector<-c("grey50","coral1","yellowgreen","yellow3","maroon3",
                                "cadetblue","brown","lightpink","steelblue4","black")

for (p in 1:10){
    empty.list[[p]]$genotype<-factor(empty.list[[p]]$genotype,
                          levels = c("WT", "5xFAD", "Cnp_ko","Cnp_5xFAD"))
    empty.list[[p]]$mean<-as.numeric(empty.list[[p]]$mean)
    empty.list[[p]]$sd<-as.numeric(empty.list[[p]]$sd)
    
    plot<- ggplot(empty.list[[p]], aes(x=genotype, y=mean, group=cluster, color=cluster)) + 
  geom_line() +
  geom_point()+
  ylab(paste0("Average center value (n=",cl$size[[p]],")"))+
    theme(text = element_text(size=20))+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                 position=position_dodge(0.05))+
    scale_color_manual(values = color_vector[[p]])
    
    print(plot)
    
    setEPS()
    postscript(file = paste0("/media/tsun/Data/Tsun/Constanze/2020_microglia_RNAseq/Visualization/DEG/adjp001_kmeans10/2020_Constanze_ADbulkRNA_adjp001k10_cluster",
                            p,"_TPMzscore.eps"))
    print(plot)
    dev.off()
}






