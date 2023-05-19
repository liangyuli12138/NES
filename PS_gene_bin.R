
# 2021-11-5
# plot for bin 
# PS 
#======================================================================================================================================================
library(Seurat)

ADDBIN <- function(dat,gene,freshold){
	library(Seurat)
	library(ggplot2)
	tmp.bin <- ifelse(dat[["RNA"]]@data[gene,]>freshold,1,0)
	dat@meta.data <- cbind(dat@meta.data,tmp.bin)
	colnames(dat@meta.data)[ncol(dat@meta.data)] <- paste0(gene,"_bin")
	return(dat)
}

rd001 <- readRDS("/public/workspace/lily/PS/Final_716/lesion1.RDS")

rd001 <- ADDBIN(rd001,gene="ANXA1",2)
rd001 <- ADDBIN(rd001,gene="FPR3",0.5)
rd001 <- ADDBIN(rd001,gene="CCL2",2)
rd001 <- ADDBIN(rd001,gene="CCR10",0.5)
pdf("/public/workspace/lily/PS/Final_716/gene_bin/RD001_gene_bin.pdf",useDingbats=F)
FeaturePlot(rd001,features=c('ANXA1_bin',"FPR3_bin"),blend=T,combine=F,order=T,pt.size=1,cols=c("red","blue"))[[3]]
FeaturePlot(rd001,features=c('CCL2_bin',"CCR10_bin"),blend=T,combine=F,order=T,pt.size=1,cols=c("red","blue"))[[3]]
dev.off()
#================================================================================================================================================

rd002 <- readRDS("/public/workspace/lily/PS/Final_716/lesion2.RDS")

rd002 <- ADDBIN(rd002,gene="ANXA1",2)
rd002 <- ADDBIN(rd002,gene="FPR3",0.5)
rd002 <- ADDBIN(rd002,gene="CCL2",2)
rd002 <- ADDBIN(rd002,gene="CCR10",0.5)
pdf("/public/workspace/lily/PS/Final_716/gene_bin/RD002_gene_bin.pdf",useDingbats=F)
FeaturePlot(rd002,features=c('ANXA1_bin',"FPR3_bin"),blend=T,combine=F,order=T,pt.size=1,cols=c("red","blue"))[[3]]
FeaturePlot(rd002,features=c('CCL2_bin',"CCR10_bin"),blend=T,combine=F,order=T,pt.size=1,cols=c("red","blue"))[[3]]
dev.off()
#================================================================================================================================================

dat673   = readRDS('/public/workspace/lily/PS/tiantan/P673_jc.RDS')
dat673 <- ADDBIN(dat673,gene="ANXA1",2)
dat673 <- ADDBIN(dat673,gene="FPR3",0.5)
dat673 <- ADDBIN(dat673,gene="CCL2",2)
dat673 <- ADDBIN(dat673,gene="CCR10",0.5)
pdf("/public/workspace/lily/PS/Final_716/gene_bin/P673_gene_bin.pdf",useDingbats=F)
FeaturePlot(dat673,features=c('ANXA1_bin',"FPR3_bin"),blend=T,combine=F,order=T,pt.size=1,cols=c("red","blue"))[[3]]
FeaturePlot(dat673,features=c('CCL2_bin',"CCR10_bin"),blend=T,combine=F,order=T,pt.size=1,cols=c("red","blue"))[[3]]
dev.off()
#================================================================================================================================================

dat689   = readRDS('/public/workspace/lily/PS/tiantan/P689.RDS')
dat689 <- ADDBIN(dat689,gene="ANXA1",2)
dat689 <- ADDBIN(dat689,gene="FPR3",0.5)
dat689 <- ADDBIN(dat689,gene="CCL2",2)
dat689 <- ADDBIN(dat689,gene="CCR10",0.5)
pdf("/public/workspace/lily/PS/Final_716/gene_bin/P689_gene_bin.pdf",useDingbats=F)
FeaturePlot(dat689,features=c('ANXA1_bin',"FPR3_bin"),blend=T,combine=F,order=T,pt.size=1,cols=c("red","blue"))[[3]]
FeaturePlot(dat689,features=c('CCL2_bin',"CCR10_bin"),blend=T,combine=F,order=T,pt.size=1,cols=c("red","blue"))[[3]]
dev.off()
#===============================================================================================================================================

dat912L  = readRDS('/public/workspace/lily/PS/tiantan/P912.L.RDS')
dat912L <- ADDBIN(dat912L,gene="ANXA1",2)
dat912L <- ADDBIN(dat912L,gene="FPR3",0.5)
dat912L <- ADDBIN(dat912L,gene="CCL2",2)
dat912L <- ADDBIN(dat912L,gene="CCR10",0.5)
pdf("/public/workspace/lily/PS/Final_716/gene_bin/P912L_gene_bin.pdf",useDingbats=F)
FeaturePlot(dat912L,features=c('ANXA1_bin',"FPR3_bin"),blend=T,combine=F,order=T,pt.size=1,cols=c("red","blue"))[[3]]
FeaturePlot(dat912L,features=c('CCL2_bin',"CCR10_bin"),blend=T,combine=F,order=T,pt.size=1,cols=c("red","blue"))[[3]]
dev.off()
#===============================================================================================================================================

dat912R  = readRDS('/public/workspace/lily/PS/tiantan/P912.R.RDS')
dat912R <- ADDBIN(dat912R,gene="ANXA1",2)
dat912R <- ADDBIN(dat912R,gene="FPR3",0.5)
dat912R <- ADDBIN(dat912R,gene="CCL2",2)
dat912R <- ADDBIN(dat912R,gene="CCR10",0.5)
pdf("/public/workspace/lily/PS/Final_716/gene_bin/P912R_gene_bin.pdf",useDingbats=F)
FeaturePlot(dat912R,features=c('ANXA1_bin',"FPR3_bin"),blend=T,combine=F,order=T,pt.size=1,cols=c("red","blue"))[[3]]
FeaturePlot(dat912R,features=c('CCL2_bin',"CCR10_bin"),blend=T,combine=F,order=T,pt.size=1,cols=c("red","blue"))[[3]]
dev.off()












#====================================================================================================================================================
# 2021-12-31 
# add H3 H4 RDS
ADDBIN <- function(dat,gene,freshold){
	library(Seurat)
	library(ggplot2)
	tmp.bin <- ifelse(dat[["RNA"]]@data[gene,]>freshold,1,0)
	dat@meta.data <- cbind(dat@meta.data,tmp.bin)
	colnames(dat@meta.data)[ncol(dat@meta.data)] <- paste0(gene,"_bin")
	return(dat)
}

H3  = readRDS('/public/workspace/wulx/missions/PS/newPairs/H3.RDS')

H3 <- ADDBIN(H3,gene="ANXA1",2)
H3 <- ADDBIN(H3,gene="FPR3",0.5)
H3 <- ADDBIN(H3,gene="CCL2",2)
H3 <- ADDBIN(H3,gene="CCR10",0.5)

pdf("/public/workspace/lily/PS/H3_gene_bin.pdf",useDingbats=F)
FeaturePlot(H3,features=c('ANXA1_bin',"FPR3_bin"),blend=T,combine=F,order=T,pt.size=1,cols=c("red","blue"))[[3]]
FeaturePlot(H3,features=c('CCL2_bin',"CCR10_bin"),blend=T,combine=F,order=T,pt.size=1,cols=c("red","blue"))[[3]]
dev.off()

#=====================================================================================================================================================
# H4 
H4  = readRDS('/public/workspace/wulx/missions/PS/newPairs/H4.RDS')

H4 <- ADDBIN(H4,gene="ANXA1",2)
H4 <- ADDBIN(H4,gene="FPR3",0.5)
H4 <- ADDBIN(H4,gene="CCL2",2)
H4 <- ADDBIN(H4,gene="CCR10",0.5)

pdf("/public/workspace/lily/PS/H4_gene_bin.pdf",useDingbats=F)
FeaturePlot(H4,features=c('ANXA1_bin',"FPR3_bin"),blend=T,combine=F,order=T,pt.size=1,cols=c("red","blue"))[[3]]
FeaturePlot(H4,features=c('CCL2_bin',"CCR10_bin"),blend=T,combine=F,order=T,pt.size=1,cols=c("red","blue"))[[3]]
dev.off()





























ADDBIN <- function(dat,gene,freshold){
	library(Seurat)
	library(ggplot2)
	tmp.bin <- ifelse(dat[["RNA"]]@data[gene,]>freshold,1,0)
	if(length(grep(paste0(gene,"_bin"),colnames(dat@meta.data)))>0){
		dat@meta.data[,grep(paste0(gene,"_bin"),colnames(dat@meta.data))] <- NULL
	}
	dat@meta.data <- cbind(dat@meta.data,tmp.bin)
	colnames(dat@meta.data)[ncol(dat@meta.data)] <- paste0(gene,"_bin")
	return(dat)
}

dat <- readRDS("/public/workspace/wulx/missions/PS/inte_4.RDS")

dat <- ADDBIN(dat,gene="ANXA1",2)
dat <- ADDBIN(dat,gene="FPR3",0.5)
dat <- ADDBIN(dat,gene="CCL2",2)
dat <- ADDBIN(dat,gene="CCR10",0.5)
pdf("/public/workspace/lily/PS/Final_716/gene_bin/inte4_gene_bin.pdf",useDingbats=F)
FeaturePlot(dat,features=c('ANXA1_bin',"FPR3_bin"),blend=T,combine=F,order=T,pt.size=1,cols=c("red","blue"))[[3]]
FeaturePlot(dat,features=c('CCL2_bin',"CCR10_bin"),blend=T,combine=F,order=T,pt.size=1,cols=c("red","blue"))[[3]]
dev.off()





#==================================================================================================================
load("/public/workspace/wulx/missions/PS/multifocal_inte.RData")

inteTT001 <- ADDBIN(inteTT001,gene="ANXA1",100)
inteTT001 <- ADDBIN(inteTT001,gene="FPR3",3.16)
inteTT001 <- ADDBIN(inteTT001,gene="CCL2",100)
inteTT001 <- ADDBIN(inteTT001,gene="CCR10",3.16)
pdf("/public/workspace/lily/PS/Final_716/gene_bin/inteTT001_gene_bin.pdf",useDingbats=F)
FeaturePlot(inteTT001,features=c('ANXA1_bin',"FPR3_bin"),blend=T,combine=F,order=T,pt.size=1,cols=c("red","blue"))[[3]]
FeaturePlot(inteTT001,features=c('CCL2_bin',"CCR10_bin"),blend=T,combine=F,order=T,pt.size=1,cols=c("red","blue"))[[3]]
dev.off()


inteTT002 <- ADDBIN(inteTT002,gene="ANXA1",2000)
inteTT002 <- ADDBIN(inteTT002,gene="FPR3",500)
inteTT002 <- ADDBIN(inteTT002,gene="CCL2",2000)
inteTT002 <- ADDBIN(inteTT002,gene="CCR10",500)
pdf("/public/workspace/lily/PS/Final_716/gene_bin/inteTT002_gene_bin.pdf",useDingbats=F)
FeaturePlot(inteTT002,features=c('ANXA1_bin',"FPR3_bin"),blend=T,combine=F,order=T,pt.size=1,cols=c("red","blue"))[[3]]
FeaturePlot(inteTT002,features=c('CCL2_bin',"CCR10_bin"),blend=T,combine=F,order=T,pt.size=1,cols=c("red","blue"))[[3]]
dev.off()































