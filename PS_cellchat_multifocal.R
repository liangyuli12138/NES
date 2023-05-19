
# 2022-1-3 
# this program is used to run cellphoenDB for PS samples
#===================================================================================================================================================

#===================================================================================================================================================
P689 <- readRDS("/public/workspace/lily/PS/P689.RDS")

table(P689$seurat_clusters,P689$celltype)
table(P689$seurat_clusters,P689$type.refine)

P689$CCC.type <- "others"
P689$CCC.type[which(P689$seurat_clusters%in%c(1,3,4,6))] <- "Tumor"
P689$CCC.type[which(P689$seurat_clusters%in%c(0))] <- "Macrophage"
P689$CCC.type[which(P689$seurat_clusters%in%c(2))] <- "Oligo."
P689$CCC.type[which(P689$seurat_clusters%in%c(5))] <- "Tcell"

table(P689$type.refine,P689$CCC.type)

saveRDS(P689,file="/public/workspace/lily/PS/P689.RDS")



#=================================================================================================================================================
P673_jc <- readRDS("/public/workspace/lily/PS/P673_jc.RDS")

table(P673_jc$seurat_clusters,P673_jc$celltype)
table(P673_jc$seurat_clusters,P673_jc$type.refine)

P673_jc$CCC.type <- "others"
P673_jc$CCC.type[which(P673_jc$seurat_clusters%in%c(3,5))] <- "Tumor"
P673_jc$CCC.type[which(P673_jc$seurat_clusters%in%c(0,1,4,7,8))] <- "Macrophage"
P673_jc$CCC.type[which(P673_jc$seurat_clusters%in%c(9))] <- "Oligo."
P673_jc$CCC.type[which(P673_jc$seurat_clusters%in%c(2,6))] <- "Tcell"

table(P673_jc$type.refine,P673_jc$CCC.type)

saveRDS(P673_jc,file="/public/workspace/lily/PS/P673_jc.RDS")



#=================================================================================================================================================
P912.L <- readRDS("/public/workspace/lily/PS/P912.L.RDS")

table(P912.L$seurat_clusters,P912.L$celltype)
table(P912.L$seurat_clusters,P912.L$type.refine)

P912.L$CCC.type <- "others"
P912.L$CCC.type[which(P912.L$seurat_clusters%in%c(1,3,4,8,9))] <- "Tumor"
P912.L$CCC.type[which(P912.L$seurat_clusters%in%c(0,2,5,7))] <- "Macrophage"
P912.L$CCC.type[which(P912.L$seurat_clusters%in%c(10))] <- "Oligo."
P912.L$CCC.type[which(P912.L$seurat_clusters%in%c(6))] <- "Tcell"

table(P912.L$type.refine,P912.L$CCC.type)

saveRDS(P912.L,file="/public/workspace/lily/PS/P912.L.RDS")



#=================================================================================================================================================
P912.R <- readRDS("/public/workspace/lily/PS/P912.R.RDS")

table(P912.R$seurat_clusters,P912.R$celltype)
table(P912.R$seurat_clusters,P912.R$type.refine)

P912.R$CCC.type <- "others"
P912.R$CCC.type[which(P912.R$seurat_clusters%in%c(3,4,5))] <- "Tumor"
P912.R$CCC.type[which(P912.R$seurat_clusters%in%c(0,1,2,6,8,9))] <- "Macrophage"
P912.R$CCC.type[which(P912.R$seurat_clusters%in%c(7))] <- "Oligo."
P912.R$CCC.type[which(P912.R$seurat_clusters%in%c(6))] <- "Tcell"

table(P912.R$type.refine,P912.R$CCC.type)

saveRDS(P912.R,file="/public/workspace/lily/PS/P912.R.RDS")

#===================================================================================================================================================
#===================================================================================================================================================
#===================================================================================================================================================#===================================================================================================================================================
#===================================================================================================================================================
#===================================================================================================================================================

cellphoneDB_input <- function(mat, clusterInfo, expr_outfile, cellinfo_outfile){      
	# mat: Seurat RDS      
	# clusterInfo: group in mat@meta.data      
	count = as.data.frame(mat@assays$RNA@counts)      
	Gene = rownames(mat@assays$RNA@counts)      
	genes = as.data.frame(Gene)      
	count1 <- cbind(genes, count)      
	write.table(count1, expr_outfile,           
	row.names = FALSE, quote = FALSE, sep = "\t")      
	info <- data.frame(Cell = colnames(mat),           
	cell_type = mat@meta.data[, clusterInfo])      
	write.table(info, cellinfo_outfile,           
	row.names = FALSE, quote = FALSE, sep = "\t")  
	
}


# load samples 
library(Seurat)
for(sample in c("P689","P673_jc","P912.L","P912.R")){
	tmp <- readRDS(paste0("/public/workspace/lily/PS/",sample,".RDS"))
	tmp.sub <- subset(tmp,cells=which(tmp$CCC.type%in%c("Tumor","Macrophage","Tcell")))
	# write out 
	dir.create(paste0("/public/workspace/lily/PS/CellphoneDB/",sample,"/"))
	cellphoneDB_input(tmp.sub,'CCC.type',
		paste0("/public/workspace/lily/PS/CellphoneDB/",sample,"/",sample,"_expr.txt"),
		paste0("/public/workspace/lily/PS/CellphoneDB/",sample,"/",sample,"_cellinfo.txt")
		)

}



# and for H3 H4 samples 
H3  = readRDS('/public/workspace/wulx/missions/PS/newPairs/H3.RDS')
H4  = readRDS('/public/workspace/wulx/missions/PS/newPairs/H4.RDS')

for(sample in c("H3","H4")){
	tmp <- readRDS(paste0("/public/workspace/wulx/missions/PS/newPairs/",sample,".RDS"))
	tmp$CCC.type <- tmp$celltype4
	tmp$CCC.type[which(tmp$CCC.type=="malignant")] <- "Tumor"
	tmp$CCC.type[which(tmp$CCC.type=="myeloid")] <- "Macrophage"
	tmp$CCC.type[which(tmp$CCC.type=="oligodendrocyte")] <- "Oligo."
	tmp$CCC.type[which(tmp$CCC.type=="Tcells")] <- "Tcell"

	tmp.sub <- subset(tmp,cells=which(tmp$CCC.type%in%c("Tumor","Macrophage","Tcell")))
	# write out 
	dir.create(paste0("/public/workspace/lily/PS/CellphoneDB/",sample,"/"))
	cellphoneDB_input(tmp.sub,'CCC.type',
		paste0("/public/workspace/lily/PS/CellphoneDB/",sample,"/",sample,"_expr.txt"),
		paste0("/public/workspace/lily/PS/CellphoneDB/",sample,"/",sample,"_cellinfo.txt")
		)

}
















#=================================================================================================================================================
# and now run CellphoneDB


#====================================run in bash 
bytlib load python-3.6.6
bytlib load cellphonedb-2.1.1
#change work path
mkdir -p /public/workspace/lily/PS/CellphoneDB/P673_jc/res

cellphonedb method statistical_analysis /public/workspace/lily/PS/CellphoneDB/P673_jc/P673_jc_cellinfo.txt \
/public/workspace/lily/PS/CellphoneDB/P673_jc/P673_jc_expr.txt --threads 6 \
--output-path=/public/workspace/lily/PS/CellphoneDB/P673_jc/res --counts-data gene_name


# for P689
bytlib load python-3.6.6
bytlib load cellphonedb-2.1.1
#change work path
mkdir -p /public/workspace/lily/PS/CellphoneDB/P689/res

cellphonedb method statistical_analysis /public/workspace/lily/PS/CellphoneDB/P689/P689_cellinfo.txt \
/public/workspace/lily/PS/CellphoneDB/P689/P689_expr.txt --threads 6 \
--output-path=/public/workspace/lily/PS/CellphoneDB/P689/res --counts-data gene_name





# for P912.L
bytlib load python-3.6.6
bytlib load cellphonedb-2.1.1
#change work path
mkdir -p /public/workspace/lily/PS/CellphoneDB/P912.L/res

cellphonedb method statistical_analysis /public/workspace/lily/PS/CellphoneDB/P912.L/P912.L_cellinfo.txt \
/public/workspace/lily/PS/CellphoneDB/P912.L/P912.L_expr.txt --threads 6 \
--output-path=/public/workspace/lily/PS/CellphoneDB/P912.L/res --counts-data gene_name





# for P912.R
bytlib load python-3.6.6
bytlib load cellphonedb-2.1.1
#change work path
mkdir -p /public/workspace/lily/PS/CellphoneDB/P912.R/res

cellphonedb method statistical_analysis /public/workspace/lily/PS/CellphoneDB/P912.R/P912.R_cellinfo.txt \
/public/workspace/lily/PS/CellphoneDB/P912.R/P912.R_expr.txt --threads 6 \
--output-path=/public/workspace/lily/PS/CellphoneDB/P912.R/res --counts-data gene_name





# for H3
bytlib load python-3.6.6
bytlib load cellphonedb-2.1.1
#change work path
mkdir -p /public/workspace/lily/PS/CellphoneDB/H3/res

cellphonedb method statistical_analysis /public/workspace/lily/PS/CellphoneDB/H3/H3_cellinfo.txt \
/public/workspace/lily/PS/CellphoneDB/H3/H3_expr.txt --threads 6 \
--output-path=/public/workspace/lily/PS/CellphoneDB/H3/res --counts-data gene_name




# for H4
bytlib load python-3.6.6
bytlib load cellphonedb-2.1.1
#change work path
mkdir -p /public/workspace/lily/PS/CellphoneDB/H4/res

cellphonedb method statistical_analysis /public/workspace/lily/PS/CellphoneDB/H4/H4_cellinfo.txt \
/public/workspace/lily/PS/CellphoneDB/H4/H4_expr.txt --threads 6 \
--output-path=/public/workspace/lily/PS/CellphoneDB/H4/res --counts-data gene_name








#===================================================================================================================================================
# plot result 
# use cellchat 
library(CellChat)

for(i in c("P689","P673_jc","P912.L","P912.R","H3","H4")){

	tmp <-read.table(paste0("/public/workspace/lily/PS/CellphoneDB/",i,"/res/significant_means.txt"),header=T,sep="\t")
	count <-  matrix(apply(tmp[,13:21],2,function(x){length(which(!is.na(x)))}),nrow=3,byrow=T)
	colnames(count) <- c("Macrophage","Tcell","Tumor")
	rownames(count) <- c("Macrophage","Tcell","Tumor")
	count[1,1] <- 0
	count[2,2] <- 0
	count[3,3] <- 0

	source("~/software/cellchat_vis.R")

	pdf(paste0("/public/workspace/lily/PS/CellphoneDB/",i,"_cellchat.pdf"),useDingbats=F)
	vis(count, weight.scale = T,label.edge=T)
	dev.off()
}








# 2022-4-11
# need matrix data 
#============================================================================================================================================

library(CellChat)

cellphonedb.res <- list()
samplelist <- c("P689","P673_jc","P912.L","P912.R","H3","H4")
for(i in 1:length(samplelist)){

	tmp <-read.table(paste0("/public/workspace/lily/PS/CellphoneDB/",samplelist[i],"/res/significant_means.txt"),header=T,sep="\t")
	count <-  matrix(apply(tmp[,13:21],2,function(x){length(which(!is.na(x)))}),nrow=3,byrow=T)
	colnames(count) <- c("Macrophage","Tcell","Tumor")
	rownames(count) <- c("Macrophage","Tcell","Tumor")
	count[1,1] <- 0
	count[2,2] <- 0
	count[3,3] <- 0

	cellphonedb.res[[i]] <- count
	names(cellphonedb.res)[i] <- samplelist[i]
}














