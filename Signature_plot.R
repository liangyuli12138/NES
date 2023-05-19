

# 2022-4-7
# plot gene 
library(Seurat)
genes <- c('S100A10','FOSL2','SPP1','CAV1','ANXA1','VIM','CD44','SERPINH1','LGALS3','CEBPB','ATF5','LGALS1')

# load data 
P912.R <- readRDS("/public/workspace/lily/PS/P912.R.RDS")
P912.L <- readRDS("/public/workspace/lily/PS/P912.L.RDS")
P673_jc <- readRDS("/public/workspace/lily/PS/P673_jc.RDS")
P689 <- readRDS("/public/workspace/lily/PS/P689.RDS")
rd001 <- readRDS("/public/workspace/lily/PS/hbrd001.rds")
rd001$CCC.type <-  car::recode(rd001$llymarker,
	" 'Tumor cell'='Tumor';
	'Myeloid'='Macrophage';
	'Oligodendrocyte'='Oligo.';
	'Endothelial'='others';
	'Fibroblast/Vascular'='others' "
	)
rd002 <- readRDS("/public/workspace/lily/PS/hbrd002.rds")
rd002$CCC.type <-  car::recode(rd002$llymarker,
	" 'Tumor cell'='Tumor';
	'Myeloid'='Macrophage';
	'Oligodendrocyte'='Oligo.';
	'Endothelial'='others';
	'T_cell'='Tcell';
	'Fibroblast/Vascular'='others' "
	)
H3  = readRDS('/public/workspace/wulx/missions/PS/newPairs/H3.RDS')
H3$CCC.type <-  car::recode(H3$celltype4,
	" 'malignant'='Tumor';
	'myeloid'='Macrophage';
	'oligodendrocyte'='Oligo.';
	'Endothelial'='others';
	'Tcells'='Tcell';
	c('fibroblast_vascular','unknown')='others' "
	)
H4  = readRDS('/public/workspace/wulx/missions/PS/newPairs/H4.RDS')
H4$CCC.type <-  car::recode(H4$celltype4,
	" 'malignant'='Tumor';
	'myeloid'='Macrophage';
	'oligodendrocyte'='Oligo.';
	'Endothelial'='others';
	'Tcells'='Tcell';
	c('fibroblast_vascular','unknown')='others' "
	)







#=======================================================================================================
# plot result 
samplelist <- c("H3" ,"H4" ,"P673_jc","P689","P912.L","P912.R","rd001","rd002")
cols <- c('#ff0000','#fbb034','#00a4e4','#89ba16','#a51890')
names(cols) <- c("Tumor","Macrophage","Oligo.","Tcell","others")

datlist <- list()
for(i in 1:length(samplelist)){
	datlist[[i]] <- get(samplelist[i])
	names(datlist)[i] <- samplelist[i]
}


percent_feature <- function(dat,genelist,group){
    res.list <- c()
    for(i in 1:length(genelist)){
        dat$tmp_gene <- ifelse(dat[["RNA"]]@data[genelist[i],]>0,"Y","N")
        if(all(dat$tmp_gene=="N")){
           res.list[[i]] <- rep(0,length=length(table(dat@meta.data[,group])))
        }else{
           res.list[[i]] <- apply(table(dat$tmp_gene,dat@meta.data[,group]),1,function(x){x/sum(x)})[,2] # change ways
        }
        
        names(res.list)[i] <- genelist[i]
    }
    return(res.list)
}


list2mat <- function(res.list){
    res.dat <- matrix(unlist(res.list),ncol=length(res.list[[1]]),byrow=T)
    rownames(res.dat) <- names(res.list)
    colnames(res.dat) <- names(res.list[[1]])
    return(res.dat)
}


for(i in 1:length(datlist)){
	# 1. plot signature gene t-sne 
	pdf(paste0("/public/workspace/lily/PS/Response/Signature/",names(datlist)[i],".pdf"),useDingbats=F)
	tmp.dat <- datlist[[i]]

	p0 <- DimPlot(tmp.dat,group.by="CCC.type",cols=cols)
	print(p0)
	DefaultAssay(tmp.dat) <- "RNA"
	for(g in 1:length(genes)){

		tmp.p <- FeaturePlot(tmp.dat,features=genes[g])
		print(tmp.p)
	}

	


	# 2. plot signature gene percentage 

	tmp.dat.list <- percent_feature(dat=tmp.dat,genelist=genes,group="CCC.type")
	tmp.dat.mat <- list2mat(tmp.dat.list)

	p_e <- pheatmap::pheatmap(tmp.dat.mat,display_numbers=T,number_color="black")
	print(p_e)


	dev.off()


}
























