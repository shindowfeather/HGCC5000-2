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


load("intergrated.all.rdata")

#combined_dataset_core
new_cluster<-read.table("cell.annotation.all.public.glio.txt",header=T,sep="\t")
#new_cluster2<-data.frame(new_cluster[,2])
#rownames(new_cluster2)<-new_cluster[,1]
integratation.combined<-AddMetaData(integratation.combined,metadata=new_cluster)

cellChat <- createCellChat(object = integratation.combined, group.by = "public.CellAnnotation",assay ="RNA")

CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
cellChat@DB <- CellChatDB.use
cellChat <- subsetData(cellChat)
future::plan("multisession", workers = 8) 
cellChat <- identifyOverExpressedGenes(cellChat)
cellChat <- identifyOverExpressedInteractions(cellChat)
cellChat <- computeCommunProb(cellChat, type = "triMean")
cellChat <- computeCommunProbPathway(cellChat)
cellChat <- aggregateNet(cellChat)

save(cellChat,file="cellchat.all.rdata")

pdf("all_cell_cell_communication.network.target.pdf")
netVisual_bubble(cellChat, sources.use = c(1,2,4,11:15),targets.use = c(5:10),remove.isolate = FALSE)
netVisual_bubble(cellChat, sources.use = c(1,2,4,11:15),targets.use = c(5:10),remove.isolate = FALSE,signaling = c("OSM","MIF","SPP1"))
dev.off()

library(NMF)
pdf("nmf.all.outgoing.incoming.pdf")
selectK(cellChat, pattern = "incoming")
selectK(cellChat, pattern = "outgoing")
dev.off()

nPatterns = 3
cellChat<-identifyCommunicationPatterns(cellChat,pattern = "incoming", k = nPatterns)
cellChat<-identifyCommunicationPatterns(cellChat,pattern = "outgoing", k = nPatterns)
pdf("all.dot.communication.plot.pdf")
netAnalysis_dot(cellChat, pattern = "incoming")
netAnalysis_dot(cellChat, pattern = "outgoing")
dev.off()

pathways.show<-"MIF"
cellChat <- netAnalysis_computeCentrality(cellChat, slot.name = "netP") 
pdf("all.mif.pathway.single.pdf")
netAnalysis_signalingRole_network(cellChat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
dev.off()

nPatterns = 3
colorvector<-c("#941B2A","#AF4137","#A37660","#1b5077","#4887ad","#8b97bb")
color.use <- setNames(colorvector,c("core1","core2","core3","edge1","edge2","edge3"))
pdf("all.dot.communication.plot.selected.pdf")
netAnalysis_dot(cellChat, pattern = "incoming",group.show = c("core1","core2","core3","edge1","edge2","edge3"),color.use=color.use)
dev.off()
