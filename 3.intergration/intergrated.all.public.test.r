
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
set.seed(1234)
#library(dplyr)
library(AUCell)
library(GSEABase)
library(ggplot2)
library(future)
library(tidyverse)
library(RCurl)
library(AnnotationHub)
library(ggplot2)
library(ArchR)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(ensembldb)

#library(BiocParallel)
#plan("multicore", workers = 16)

data2<-get(load("glio.all.merged_subclusters.new.rdata"))
load("publid2.rdata")
subset_seurat<-data2

DefaultAssay(subset_seurat) <- "CCA_RNA"
features <- SelectIntegrationFeatures(object.list = list(subset_seurat,public),assay=c("CCA_RNA","RNA"))
integratation <- FindIntegrationAnchors(object.list =list(subset_seurat,public) ,scale=F, anchor.features = features,assay=c("CCA_RNA","RNA"),dims=1:30,k.anchor=500,k.score = 500)
integratation.combined <- IntegrateData(anchorset = integratation,new.assay.name = "integrated",normalization.method = c("LogNormalize"),k.weight = 500)
#integratation.combined <- RunUMAP(integratation.combined, reduction = "integrated", dims = 1:30)

DefaultAssay(integratation.combined) <- "integrated"
integratation.combined <- ScaleData(integratation.combined, verbose = FALSE)
integratation.combined <- RunPCA(integratation.combined, npcs = 30, verbose = FALSE)
integratation.combined <- RunUMAP(integratation.combined, reduction = "pca", dims = 1:30)


sample_name<-data.frame(gsub("10_GBM_616_|11_GBM_616_|12_GBM_616_|9_GBM_604_|617_|14_GBM_617_|2_GBM_609_|3_GBM_609_|4_GBM_622_|5_GBM_622_|6_GBM_622_|7_GBM_604_|13_GBM_|8_GBM_604_|15_GBM_617_","",integratation.combined$Sample))
 colnames(sample_name)<-"new_sample_name"
integratation.combined<-AddMetaData(integratation.combined,metadata=sample_name)

pdf("merged.all_merged_clusters.immune.cell.pdf")
DimPlot(integratation.combined, reduction = "umap", group.by = "CellAnnotation",pt.size=0.3)
DimPlot(integratation.combined, reduction = "umap", group.by = "sample",pt.size=0.3)
DimPlot(integratation.combined, reduction = "umap", group.by = "new_sample_name",pt.size=0.3)
dev.off()

save(integratation.combined,file="intergrated.all2.rdata")
#write.table(data.frame(integratation.combined$core_edge_cluster),file="cell.annotation.glio.txt",quote=F,sep="\t")

colorvector<-c("#941B2A","#AF4137","#A37660","#1b5077","#4887ad","#8b97bb","#d9d9d9","#bdbdbd","#969696","#737373","#6b6969","#525050","#252525","#2b2b2b","#323332")
new_cluster<-read.table("cell.annotation.all.public.glio.txt",header=T,sep="\t")
integratation.combined<-AddMetaData(integratation.combined,metadata=new_cluster)
integratation.combined$public.CellAnnotation<-factor(integratation.combined$public.CellAnnotation,levels=unique(integratation.combined$public.CellAnnotation))
pdf("merged.all_merged_clusters.immune.cell.annotation.pdf")
DimPlot(integratation.combined, reduction = "umap", group.by = "public.CellAnnotation",pt.size=0.3,cols=colorvector)
dev.off()


colorvector<-c("#ffd3e7","#ffd3e7","#ffd3e7","#adcde4","#adcde4","#adcde4","#d9d9d9","#bdbdbd","#969696","#737373","#6b6969","#525050","#252525","#2b2b2b","#323332","#6c6c6c","#949494")
colorvector<-c("#525050","#d9d9d9","#252525","#2b2b2b","#bdbdbd","#6b6969","#969696","#737373","#323332","#6c6c6c","#949494","#adcde4","#adcde4","#adcde4","#ffd3e7","#ffd3e7","#ffd3e7")
new_cluster<-read.table("cell.annotation.all.public.glio.cd45.txt",header=T,sep="\t")
integratation.combined<-AddMetaData(integratation.combined,metadata=new_cluster)
integratation.combined$public.CellAnnotation<-factor(integratation.combined$public.CellAnnotation,levels=unique(integratation.combined$public.CellAnnotation))
pdf("merged.all_merged_clusters.immune.cell.annotation.pdf")
DimPlot(integratation.combined, reduction = "umap", group.by = "public.CellAnnotation",pt.size=0.3,cols=colorvector,order="public.CellAnnotation")
dev.off()
