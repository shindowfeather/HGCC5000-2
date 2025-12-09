library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
set.seed(1234)
library(dplyr)
library(AUCell)
library(GSEABase)
library(ggplot2)
library(gridExtra)
library(viridis)
library(MASS)
library(ggpubr)
library(doRNG)

load("/proj/sens2023512/nobackup/wharf/lucycy/lucycy-sens2023512/database/list/list_mes.immune.list.rdata")

a<-get(load("/proj/sens2023512/nobackup/wharf/lucycy/lucycy-sens2023512/glioblastoma_single_cell/1.RNA/1.RNA_seq_seurat/glio.int.all.peaks.chromvar.rdata"))

core_matrix = as.matrix(a@assays$RNA@data)
core_rankings <- AUCell_buildRankings(core_matrix, nCores=30, plotStats=FALSE)
rm(core_matrix)

list2<-list()
for (i in 1:7){
     list2[[i]]<- AUCell_calcAUC(list_mes[[i]], core_rankings)
}

for(each in 1:7){
    AUC <- t(as.data.frame(list2[[each]]@assays@data$AUC))
    colnames(AUC) <- paste0(colnames(AUC), "_AUC")
    AUC <- data.frame(AUC)
    #print(AUC)
    a <- AddMetaData(a, metadata = AUC)
}

g<-a@meta.data[,c(41:47)]

annotation<-data.frame(a@meta.data$`CCA_RNA_snn_res.0.1`)
rownames(annotation)<-rownames(a@meta.data)
dim(annotation)
dim(g)

annotation2<-data.frame(annotation[order(annotation$`a.meta.data.CCA_RNA_snn_res.0.1`),])
rownames(annotation2)<-rownames(annotation)[order(annotation$`a.meta.data.CCA_RNA_snn_res.0.1`)]
colnames(annotation2)<-"group"

library(pheatmap)
library( "RColorBrewer" )
color = colorRampPalette(brewer.pal(3,"RdBu"))(21)
breaks<-seq(-2,2,by=0.15)
g2<-t(g)[,order(annotation$`a.meta.data.CCA_RNA_snn_res.0.1`)]

pheatmap(g2,show_rownames=T, annotation=annotation2,cluster_rows=F,cluster_cols = F,clustering_distance_cols="euclidean",scale="row",show_colnames= F,color = rev(color),clustering_distance_rows="euclidean", cex=1, clustering_method="complete", border_color=FALSE,fontsize_row=12,fontsize_col=12,breaks=breaks)


d<-data.frame(scale(g))
d$CCA_RNA_snn_res.0.1<-a$CCA_RNA_snn_res.0.1
d$type<-a$type
bar3<-data.frame(value=d$Hara_OSM_AUC,name=d$`CCA_RNA_snn_res.0.1`)
my_comparisons<-list(c("0","1"))
ggplot(bar3,aes(x=name,y=value,fill=name))+geom_boxplot()+scale_fill_manual(values = c("#941B2A","#AF4137","#A37660"))+theme(text=element_text(size=20,family="Helvetica"))+ylab("Z-score")+theme_minimal()+theme(text = element_text(size = 20),panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1))+stat_compare_means(comparisons = my_comparisons,label="p.signif",method="wilcox.test",paired=F)

pdf("core.edge.mes.pdf")
for (i in 1:7){
    bar3<-data.frame(value=d[,i],name=d$type)
my_comparisons<-list(c("core","edge"))
        aa<-ggplot(bar3,aes(x=name,y=value,fill=name))+geom_boxplot()+scale_fill_manual(values = c("#D4CB92","#395C6B"))+theme(text=element_text(size=20,family="Helvetica"))+ylab("Z-score")+theme_classic()+
    ggtitle(colnames(d)[i])+
    theme(text = element_text(size = 20),panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1))+
	stat_compare_means(comparisons = my_comparisons,label="p.format",method="wilcox.test",paired=F)
    print (aa)
}
dev.off()
