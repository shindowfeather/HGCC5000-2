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

load("gbm.annotation.rdata")
load("/proj/sens2023512/nobackup/wharf/lucycy/lucycy-sens2023512/database/list/single/list.single.list.rdata")

cell_type<-read.table("final1.annotation.txt",sep="\t")
cell_type2<-data.frame(cell_annotation=cell_type$V2)
rownames(cell_type2)<-cell_type$V1
gbm<-AddMetaData(gbm,metadata=cell_type2)

Idents(gbm) <- "cell_annotation"
new<-subset(gbm, idents = "tumor")

Idents(new) <- "tumor"
edge_cell<-"Peritumoral"
core_cell<-"Tumoral"
edge_tis <- subset(new, idents = edge_cell)
core_tis <- subset(new, idents = core_cell)


core_matrix = as.matrix(core_tis@assays$RNA@data)
core_rankings <- AUCell_buildRankings(core_matrix, nCores=30, plotStats=FALSE)
rm(core_matrix)

MES2_core_AUC <- AUCell_calcAUC(list$single.mes_new2, core_rankings)
MES1_core_AUC <- AUCell_calcAUC(list$single.mes_new1, core_rankings)
NPC1_core_AUC <- AUCell_calcAUC(list$single.npc_new1, core_rankings)
NPC2_core_AUC <- AUCell_calcAUC(list$single.npc_new2, core_rankings)
AC_core_AUC <- AUCell_calcAUC(list$`single.ac-1`, core_rankings)
OPC_core_AUC <- AUCell_calcAUC(list$single.opc_new, core_rankings)

core_AUC_list = c('MES1_core_AUC','MES2_core_AUC','NPC1_core_AUC',
                              'NPC2_core_AUC','AC_core_AUC','OPC_core_AUC')

for(each in core_AUC_list){
    AUC <- t(as.data.frame(get(each)@assays@data$AUC))
    colnames(AUC) <- paste0(colnames(AUC), "_AUC")
    AUC <- data.frame(AUC)
    #print(AUC)
    core_tis <- AddMetaData(core_tis, metadata = AUC)
}

### Define D (aka y-axis)
core_meta <- core_tis@meta.data

core_meta$MES_mean <- apply(cbind(core_meta$MES1_AUC,core_meta$MES2_AUC),1,mean)
core_meta$NPC_mean <- apply(cbind(core_meta$NPC1_AUC,core_meta$NPC2_AUC),1,mean)

#take max of these two scores
a1 <- apply(cbind(core_meta$OPC_AUC, core_meta$NPC_mean), 1, max)
b1 <- apply(cbind(core_meta$AC_AUC, core_meta$MES_mean), 1, max)
core_meta$Y.axis <-  a1 - b1
core_meta$CellClass <- ifelse(core_meta$Y.axis > 0, "OPC.NPC", "AC.MES")

head(core_meta)

### classify cells by their max score

cc <- cbind(core_meta$AC_AUC, core_meta$MES_mean, core_meta$NPC_mean, core_meta$OPC_AUC)
colnames(cc) <- c("AC", "MES", "NPC", "OPC")
core_meta$MaxClass <- apply(cc,1,function(x) which(x==max(x)))
core_meta$MaxClass <- gsub(1, "AC", core_meta$MaxClass)
core_meta$MaxClass <- gsub(2, "MES", core_meta$MaxClass)
core_meta$MaxClass <- gsub(3, "NPC", core_meta$MaxClass)
core_meta$MaxClass <- gsub(4, "OPC", core_meta$MaxClass)
    
core_meta$ClassSign_X <- core_meta$MaxClass
core_meta$ClassSign_X <- gsub("AC", -1, core_meta$ClassSign_X )
core_meta$ClassSign_X <- gsub("MES", 1,core_meta$ClassSign_X )
core_meta$ClassSign_X <- gsub("OPC", -1, core_meta$ClassSign_X )
core_meta$ClassSign_X <- gsub("NPC", 1, core_meta$ClassSign_X )

core_meta$OPC.NPC_x <- log2(abs(core_meta$OPC_AUC - core_meta$NPC_mean)+1)
core_meta$AC.MES_x <- log2(abs(core_meta$AC_AUC - core_meta$MES_mean)+1)
core_meta$X.axis <- NA
core_meta$X.axis[grep("OPC.NPC", core_meta$CellClass)] <- core_meta$OPC.NPC_x[grep("OPC.NPC", core_meta$CellClass)]
core_meta$X.axis[grep("AC.MES", core_meta$CellClass)] <- core_meta$AC.MES_x[grep("AC.MES", core_meta$CellClass)]
#now change the sign based on maxclass for cell
#if AC is the max value, make postive; if MES is the max, make negative
#OPC = +, NPC = -
core_meta$X.axis_Class <- core_meta$X.axis * as.numeric(core_meta$ClassSign_X)

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

core_meta$density <- get_density(core_meta$X.axis_Class, core_meta$Y.axis, n = 100)

d <- ggplot(core_meta, aes(x = X.axis_Class, y = Y.axis, color = density)) +
geom_point(aes(X.axis_Class, Y.axis, color = density), alpha = 1, size = 0.7) +  scale_color_viridis() + 
xlab("<-- AC-like      MES-like -->") +
ylab("<-- AC-like  OPC-like -->") + 
geom_vline(xintercept = 0, lty = 2, col = "black") + 
geom_hline(yintercept = 0, lty = 2, col = "black") + 
ggtitle("<-- OPC-like      NPC-like -->")  +
 theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))+ylim(-0.2,0.2)+xlim(-0.4,0.4)
d

library(RColorBrewer)
core_meta$density <- get_density(core_meta$X.axis_Class, core_meta$Y.axis, n = 100)

core_meta<-as.data.frame(core_meta)
#core_meta$X.axis_Class <- factor(core_meta$X.axis_Class)

d <- ggplot(core_meta, aes(x = X.axis_Class, y = Y.axis, color = density)) +
geom_point(aes(X.axis_Class, Y.axis, color = density), alpha = 1, size = 0.7) + geom_density_2d(bins = 3, color = "grey")+scale_colour_gradientn(colours = rev(colorRampPalette(brewer.pal(11,"Greys"))(31))) +
xlab("<-- AC-like      MES-like -->") +
ylab("<-- AC-like  OPC-like -->") + 
geom_vline(xintercept = 0, lty = 2, col = "black") + 
geom_hline(yintercept = 0, lty = 2, col = "black") + 
ggtitle("<-- OPC-like      NPC-like -->")  +
 theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))+ylim(-0.3,0.3)+xlim(-0.5,0.5)
d

a<-table(core_meta$MaxClass,core_meta$tumor )

new<-t(t(a)/colSums(a))
new2<-data.frame(new)
#new2

new2$Var1<-factor(new2$Var1,levels=c("NPC","OPC","AC","MES"))
new2$sample<-c(rep("core",4))
pdf("core.barplot.pdf")
ggplot(new2, aes( y=Freq, x=Var2))+scale_fill_manual(values=c("#0571B0","#92C5DE","#F4A582","#CA0020"))+geom_col(aes(fill = Var1 ), width = 0.7)+theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

write.table(new2,file="new.core.txt",quote=F)
