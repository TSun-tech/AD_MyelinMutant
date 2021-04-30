#2020, Ting Sun
#for Depp and Sun et al manuscript
#microglia bulk RNA sequencing
#raw fastq files aligned to reference genome using STAR default parameters
#gene counts extracted by featureCounts

#sequenced genotypes (6-month-old)
#WT, CnpKO, 5xFAD, CnpKO x 5xFAD

#apply differential expression analysis to every pair of genotypes (using CnpKOx5xFAD vs CnpKO as example )
#by DESeq2 package
#summarize statistic result


library(dplyr)
library(DESeq2)
library(ggplot2)

sessionInfo()

#read in raw counts matrix, sample meta.data list and gene annotation list
#meta.data with sample ID, path, matched to genotype and replicate number
#annot table with Ensemble ID matched to gene symbol (see counts_TPM_PCA script)

meta.data<-read.csv("./bulk_metadata.csv")
annot<-read.csv("./mouse_gene_ID_annotation.csv",
               stringsAsFactors = FALSE)
raw<-read.csv("./2020_AD_microglia_raw_counts.csv")

meta.data
head(annot)
head(raw)
colnames(raw)

rownames(raw)<-raw$X
rownames(meta.data)<-meta.data$sample

colnames(raw)<-gsub("X","",colnames(raw))
colnames(raw)



#####################select pair genotypes for analysis
#examplie using CnpKOx5xFAD vs CnpKO
sampleTable <- data.frame(sample = meta.data$sample[c(5:8,13:16)],
                         condition = meta.data$genotype[c(5:8,13:16)])
rownames(sampleTable)<-sampleTable$sample

sampleTable

colnames(raw)[1]<-"X"
rownames(raw)<-raw$X
cts<-raw[,rownames(sampleTable)]
head(cts)

#apply DGE analysis using DESeq2 package
dds<-DESeqDataSetFromMatrix(countData = cts,
                                    colData = sampleTable,
                                    design = ~condition)

dds

##############
dds$condition<-relevel(dds$condition, ref = "Cnp_ko")

dds_result<-DESeq(dds)

#extract statistic result
res<-results(dds_result)

res

exp<-res
########################
colnames(exp)<-paste0("Cnp5xFAD_vs_Cnpko_",colnames(res))

exp

#order result based on adjP
res<-res[order(res$padj, decreasing = F),]

#filter result for significant regulated genes, use cut-off adjP<0.05
res
summary(res)

res_sig<-results(dds_result, alpha = 0.05)
summary(res_sig)

plotMA(res_sig, ylim = c(-10,10))

mcols(res_sig)


#########################the following script aim to summarize statistics from CnpKOx5xFAD vs CnpKO
##########################together with gene raw counts and normalised counts

raw<-as.data.frame(counts(dds_result,normalized=F))
rld<-as.data.frame(counts(dds_result,normalized=T))

head(raw)
head(rld)

for (i in 1:8){
    colnames(raw)[[i]]<-paste0(colnames(raw)[[i]],".","raw_counts")
    colnames(rld)[[i]]<-paste0(colnames(rld)[[i]],".","_normalised_value")
}

head(raw)
head(rld)

raw$gene_name<-rownames(raw)
rld$gene_name<-rownames(rld)
res$gene_name<-rownames(res)

res<-as.data.frame(res[rownames(raw),])

################################
colnames(res)[1:6]<-paste0("Cnp5xFAD_vs_Cnpko_",colnames(res)[1:6])

pre<-merge(x = res, y = raw, by = "gene_name")
all<-merge(x = pre, y = rld, by = "gene_name")

head(all)

#######################
all<-all[order(all$Cnp5xFAD_vs_Cnpko_padj, decreasing = F),]

head(all)

head(annot)

all$ensembl_ID<-all$gene_name

for(i in 1:nrow(all)){
    position<-match(all$ensembl_ID[[i]], annot$ensembl_gene_id)
    all$gene_name[[i]]<-annot$mgi_symbol[position]
}

head(all)

#just for OCD
all<-all[,c(24,1:23)]

head(all)
