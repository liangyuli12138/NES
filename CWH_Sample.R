

# 2022-5-4
# analysis CWH three samples
# prepare 
#========================================================================================================================================
library(Seurat)
for (sample in c("CWH_1","CWH_2","CWH_3")){
	filepath = paste0("/public/workspace/lily/PS/Response/",sample,"/outs/filtered_feature_bc_matrix")
	respath = "/public/workspace/lily/PS/Response/CWH_analysis/"


	tmp <- Read10X(data.dir = filepath)
	dat <- CreateSeuratObject(counts = tmp,  project = sample,min.cells = 3, min.features = 200)
	dat[["percent.mt"]] <- PercentageFeatureSet(object = dat, pattern = "^MT-")
	pdf(paste0(respath,sample, "_vlnPlot_prepare.pdf"))
	p <- VlnPlot(object = dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
	print(p)
	dev.off()



	dat = subset(x=dat,subset=nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
	

	#========== normalize and find cluster
	tmp_dat = NormalizeData(object = dat)
	tmp_dat <- FindVariableFeatures(object = tmp_dat)
	# scaling
	all.genes <- rownames(x = tmp_dat)
	tmp_dat <- ScaleData(object = tmp_dat, features = all.genes)
	# PCA
	tmp_dat <- RunPCA(object = tmp_dat, features = VariableFeatures(object = tmp_dat))
	# clustering
	tmp_dat <- FindNeighbors(object = tmp_dat,dims=1:10)
	# select proper resolution
	tmp_dat <- FindClusters(object = tmp_dat,resolution=0.8)
	# T-SNE
	tmp_dat <- RunTSNE(object = tmp_dat,dims=1:10)
	tmp_dat <- RunUMAP(tmp_dat,dims=1:10)


	saveRDS(tmp_dat,file=paste0(respath,sample,'.RDS'))


}




# check if have tumor cells 
#=============================================================================================================================================

library(Seurat)
for (sample in c("CWH_1","CWH_2","CWH_3")){
	filepath = paste0("/public/workspace/lily/PS/Response/CWH_analysis/",sample,".RDS")
	respath = paste0("/public/workspace/lily/PS/Response/CWH_analysis/",sample,"_")
	dat <- readRDS(filepath)


	p1 = DimPlot(dat,label=T,label.size=6)
	p2 = FeaturePlot(dat,features=c('CD3D','CD3E','CD2','PTPRC'),label=T,label.size=3,order=T) # T cell
	p3 = FeaturePlot(dat,features=c('CD19','CD68','FCGR3A','LYZ'),label=T,label.size=3,order=T) # Myeloid
	p4 = FeaturePlot(dat,features=c('MS4A1',"CD79A",'PTPRC'),label=T,label.size=3,order=T) # B cell 
	p5 = FeaturePlot(dat,features=c('MAG','MOG','CNDP1','PTPRC'),label=T,label.size=3,order=T) # Oligodendrocyte
	p6 = FeaturePlot(dat,features=c('COL1A1','COL1A2','DCN','CD248'),label=T,label.size=3,order=T) # Fibroblast/Vascular
	p7 = FeaturePlot(dat,features=c('CLDN5','VWF','ABCG2','CDH5'),label=T,label.size=3,order=T) # Endothelial
	p8 = FeaturePlot(dat,features=c("EGFR","PTPRZ1","SOX2"),label=T,label.size=3,order=T,cols=c("lightgrey", "red"))



	p9 = DimPlot(dat,label=T,label.size=6)
	p10 = VlnPlot(dat,features=c('CD3D','CD3E','CD2','PTPRC'),pt.size=0) # T cell
	p11 = VlnPlot(dat,features=c('CD19','CD68','FCGR3A','LYZ'),pt.size=0) # Myeloid
	p12 = VlnPlot(dat,features=c('MS4A1',"CD79A",'PTPRC'),pt.size=0) # B cell 
	p13 = VlnPlot(dat,features=c('MAG','MOG','CNDP1','PTPRC'),pt.size=0) # Oligodendrocyte
	p14 = VlnPlot(dat,features=c('COL1A1','COL1A2','DCN','CD248'),pt.size=0) # Fibroblast/Vascular
	p15 = VlnPlot(dat,features=c('CLDN5','VWF','ABCG2','CDH5'),pt.size=0) # Endothelial
	p16 = VlnPlot(dat,features=c("EGFR","PTPRZ1","SOX2"),pt.size=0)



	pdf(paste0(respath,"FeaturePlot.pdf"),useDingbats=F)
	print(p1)
	print(p2)
	print(p3)
	print(p4)
	print(p5)
	print(p6)
	print(p7)
	print(p8)
	dev.off()


	pdf(paste0(respath,"VlnPlot.pdf"),useDingbats=F,width=12)
	print(p9)
	print(p10)
	print(p11)
	print(p12)
	print(p13)
	print(p14)
	print(p15)
	print(p16)
	dev.off()







}







#===============================================================================================================================================
# celltype analysis 
# 2022-5-4

CWH1 <- readRDS("/public/workspace/lily/PS/Response/CWH_analysis/CWH_1.RDS")
CWH2 <- readRDS("/public/workspace/lily/PS/Response/CWH_analysis/CWH_2.RDS")
CWH3 <- readRDS("/public/workspace/lily/PS/Response/CWH_analysis/CWH_3.RDS")


CWH1$celltype <-  car::recode(CWH1$seurat_clusters,
	" c('16','4','20')='Tcell';
	c('1','2','11','12','13','15','17','21','6')='Macrophage';
	c('8','18')='Endothelial';
	'19'='fibroblast_vascular';
	c('0','3','5','7','9','10','14')='Tumor';
	"
	)


CWH2$celltype <-  car::recode(CWH2$seurat_clusters,
	" c('9','5')='Tcell';
	c('0','1','2','3','4','8','13','16')='Macrophage';
	c('6','15')='Endothelial';
	'17'='Oligo.';
	'12'='fibroblast_vascular';
	c('7','11','10','14')='Tumor';
	"
	)

CWH3$celltype <-  car::recode(CWH3$seurat_clusters,
	" c('12')='Tcell';
	c('9','13','14')='Macrophage';
	'15'='fibroblast_vascular';
	c('0','1','2','3','4','5','6','7','8','10','11')='Tumor';
	"
	)



saveRDS(CWH1,file="/public/workspace/lily/PS/Response/CWH_analysis/CWH_1.RDS")
saveRDS(CWH2,file="/public/workspace/lily/PS/Response/CWH_analysis/CWH_2.RDS")
saveRDS(CWH3,file="/public/workspace/lily/PS/Response/CWH_analysis/CWH_3.RDS")









#==============================================================================================================================================
# 2022-5-4
# run infercnv
library(infercnv)

CWH1 <- readRDS("/public/workspace/lily/PS/Response/CWH_analysis/CWH_1.RDS")
CWH2 <- readRDS("/public/workspace/lily/PS/Response/CWH_analysis/CWH_2.RDS")
CWH3 <- readRDS("/public/workspace/lily/PS/Response/CWH_analysis/CWH_3.RDS")




	# need to change 
	i = "CWH_3"
	CWH3$CNV.type <- car::recode(CWH3$celltype,
	" c('Tcell','Macrophage','Endothelial','Oligo.','fibroblast_vascular')='nonTumor';
	c('Tumor')='Tumor';
	"
	)

    tmp.dat <- CWH3

    DefaultAssay(tmp.dat) <- "RNA"
    dir.create(paste0("/public/workspace/lily/PS/Response/CWH_analysis/InferCNV/",i))
    setwd(paste0("/public/workspace/lily/PS/Response/CWH_analysis/InferCNV/",i))
    # respath <- paste0("/public/workspace/lily/Lung2Brain/Version6/Prepare/InferCNV/",samplelist[i])

    write.table(tmp.dat$CNV.type,paste(i,"_cell_info.txt",sep='_'),sep="\t",col.names=F,quote=F)
    count <- as.matrix(tmp.dat@assays$RNA@counts)
    write.table(count,paste(i,"_count_exp.txt",sep='_'),sep="\t",quote=F)

    infercnv_obj = CreateInfercnvObject(raw_counts_matrix=paste(i,"_count_exp.txt",sep='_'),
            annotations_file=paste(i,"_cell_info.txt",sep='_'),
            delim="\t",
            gene_order_file="/public/workspace/lily/REF/INDEX-hg19/anno/gencode_hg19_pos.txt",
            ref_group_names=c("nonTumor"))
            
    infercnv_obj = infercnv::run(infercnv_obj,
                                cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                out_dir="./", 
                                cluster_by_groups=T, 
                                plot_steps=F,
                                no_prelim_plot = TRUE,
                                num_threads=8, #big
                                no_plot=F ,
                                output_format = "pdf" # maybe can more quick 
                                # used for final scaling to fit range (0,2) centered at 1.
                                )













#======================================================================================================================================================
# 2022-5-25
# run velocity for these samples 
# run 202.195.187.3

conda activate velocity
module load  samtools-1.9

velocyto run10x -m /public/workspace/lily/REF/hg19_rmsk.gtf \
    -@ 10 --samtools-memory 2000 \
    /public/workspace/lily/PS/Response/CWH_1/ \
    /public/workspace/lily/REF/refdata-cellranger-hg19-1.2.0/genes/genes.gtf




conda activate velocity
module load  samtools-1.9

velocyto run10x -m /public/workspace/lily/REF/hg19_rmsk.gtf \
    -@ 10 --samtools-memory 2000 \
    /public/workspace/lily/PS/Response/CWH_2/ \
    /public/workspace/lily/REF/refdata-cellranger-hg19-1.2.0/genes/genes.gtf




conda activate velocity
module load  samtools-1.9

velocyto run10x -m /public/workspace/lily/REF/hg19_rmsk.gtf \
    -@ 10 --samtools-memory 2000 \
    /public/workspace/lily/PS/Response/CWH_3/ \
    /public/workspace/lily/REF/refdata-cellranger-hg19-1.2.0/genes/genes.gtf












#====================================================================================================================================================
# run result 
# Seurat object use to tumor loom 
bytlib load hdf5-1.8.13
bytlib load R-3.6.0
library(Seurat)
library(loomR)
library(velocyto.R)
library(pagoda2)
library(SCopeLoomR)
library(SeuratWrappers)
library(rlist)


loomfile = "/public/workspace/lily/PS/Response/CWH/loom_file/CWH_1.loom"
RDSfile = "/public/workspace/lily/PS/Response/CWH/CWH_1.RDS"
SampleName = "CWH1"

loomdata <- as.Seurat(ReadVelocity(loomfile))
loomdata$cellname <- gsub("x$|^.*:","",colnames(loomdata))
dat.obj <- readRDS(RDSfile)
dat.obj$cellname <- gsub("\\.1$","",colnames(dat.obj))
tumor.dat <- subset(dat.obj,cells=which(dat.obj$celltype=="Tumor"))


# subset cell 
loomdata.subset <- subset(loomdata,cells=which(loomdata$cellname%in%tumor.dat$cellname))
loomdata.subset = NormalizeData(loomdata.subset, verbose = FALSE)
loomdata.subset = FindVariableFeatures(loomdata.subset, selection.method = "vst", verbose = FALSE)
loomdata.subset <- RenameCells(loomdata.subset,new.names= paste0(SampleName,loomdata.subset$cellname))

# some tumor cell do not found ,try to subset seurat obj 
ncol(loomdata.subset) == length(dat.obj$celltype=="Tumor")
saveRDS(loomdata.subset,file="/public/workspace/lily/PS/Response/CWH/loom_file/CWH1.tumor.loom.RDS")


# for CWH2
#=====================================================================================================================
loomfile = "/public/workspace/lily/PS/Response/CWH/loom_file/CWH_2.loom"
RDSfile = "/public/workspace/lily/PS/Response/CWH/CWH_2.RDS"
SampleName = "CWH2"

loomdata <- as.Seurat(ReadVelocity(loomfile))
loomdata$cellname <- gsub("x$|^.*:","",colnames(loomdata))
dat.obj <- readRDS(RDSfile)
dat.obj$cellname <- gsub("\\.1$","",colnames(dat.obj))
tumor.dat <- subset(dat.obj,cells=which(dat.obj$celltype=="Tumor"))


# subset cell 
loomdata.subset <- subset(loomdata,cells=which(loomdata$cellname%in%tumor.dat$cellname))
loomdata.subset = NormalizeData(loomdata.subset, verbose = FALSE)
loomdata.subset = FindVariableFeatures(loomdata.subset, selection.method = "vst", verbose = FALSE)
loomdata.subset <- RenameCells(loomdata.subset,new.names= paste0(SampleName,loomdata.subset$cellname))

# some tumor cell do not found ,try to subset seurat obj 
ncol(loomdata.subset) == length(which(dat.obj$celltype=="Tumor"))
saveRDS(loomdata.subset,file=paste0("/public/workspace/lily/PS/Response/CWH/loom_file/",SampleName,".tumor.loom.RDS"))




# for CWH3
#=====================================================================================================================
loomfile = "/public/workspace/lily/PS/Response/CWH/loom_file/CWH_3.loom"
RDSfile = "/public/workspace/lily/PS/Response/CWH/CWH_3.RDS"
SampleName = "CWH3"

loomdata <- as.Seurat(ReadVelocity(loomfile))
loomdata$cellname <- gsub("x$|^.*:","",colnames(loomdata))
dat.obj <- readRDS(RDSfile)
dat.obj$cellname <- gsub("\\.1$","",colnames(dat.obj))
tumor.dat <- subset(dat.obj,cells=which(dat.obj$celltype=="Tumor"))


# subset cell 
loomdata.subset <- subset(loomdata,cells=which(loomdata$cellname%in%tumor.dat$cellname))
loomdata.subset = NormalizeData(loomdata.subset, verbose = FALSE)
loomdata.subset = FindVariableFeatures(loomdata.subset, selection.method = "vst", verbose = FALSE)
loomdata.subset <- RenameCells(loomdata.subset,new.names= paste0(SampleName,loomdata.subset$cellname))

# some tumor cell do not found ,try to subset seurat obj 
ncol(loomdata.subset) == length(which(dat.obj$celltype=="Tumor"))
saveRDS(loomdata.subset,file=paste0("/public/workspace/lily/PS/Response/CWH/loom_file/",SampleName,".tumor.loom.RDS"))





#========================================================================================================================================
CWH1.loomdata.subset <- readRDS("/public/workspace/lily/PS/Response/CWH/loom_file/CWH1.tumor.loom.RDS")
CWH2.loomdata.subset <- readRDS("/public/workspace/lily/PS/Response/CWH/loom_file/CWH2.tumor.loom.RDS")

loomDat <- merge(CWH1.loomdata.subset,CWH2.loomdata.subset)
saveRDS(loomDat,file='/public/workspace/lily/PS/Response/CWH/loom_file/CWH12.tumor.loom.rds')


# and now integration tumor seurat obj to get umap info 
dat.obj.CWH1 <- readRDS("/public/workspace/lily/PS/Response/CWH/CWH_1.RDS")
dat.obj.CWH1$cellname <- gsub("\\.1$","",colnames(dat.obj.CWH1))
tumor.dat.CWH1 <- subset(dat.obj.CWH1,cells=which(dat.obj.CWH1$celltype=="Tumor"))

dat.obj <- readRDS("/public/workspace/lily/PS/Response/CWH/CWH_2.RDS")
dat.obj$cellname <- gsub("\\.1$","",colnames(dat.obj))
tumor.dat.CWH2 <- subset(dat.obj,cells=which(dat.obj$celltype=="Tumor"))

integration.anchors <- FindIntegrationAnchors(object.list = c(tumor.dat.CWH1,tumor.dat.CWH2))
inte <- IntegrateData(anchorset = integration.anchors)
#FindVariableFeatures
inte <- FindVariableFeatures(inte)
##Scaling the integrateda
all.genes <- rownames(inte)
inte <- ScaleData(inte, features = all.genes)
#PCA
inte <- RunPCA(inte)
#cluster
inte <- FindNeighbors(inte)
inte <- FindClusters(inte)
#TSNE
# if Umap can not use
inte <- RunTSNE(inte)

inte <- RunUMAP(inte,dims=1:10)
saveRDS(inte,file="/public/workspace/lily/PS/Response/CWH/CWH12.tumor.dat.RDS")












#==================================================================================================================================================
# 2022-6-18
# run velocyto.R 

bytlib load libraries/hdf5-1.8.13
bytlib load R-3.6.0
R

library(loomR)
library(velocyto.R)
library(pagoda2)
library(SCopeLoomR)
library(SeuratWrappers)
library(rlist)
library(Seurat)




savepath="/public/workspace/lily/PS/Multifocal/Velocyto_res/CWH12/"
dat_rdsPath ="/public/workspace/lily/PS/Response/CWH/CWH12.tumor.dat.RDS"
loomDat_rdsPath ="/public/workspace/lily/PS/Response/CWH/loom_file/CWH12.tumor.loom.rds"
spliced_minAvg=0.5
unspliced_minAvg=0.1
kCells=500
n=5000
fit_quantile=0.02
param="orig.ident"

get_color_scheme = function(type = "clusters") {
  library(ggsci)
  if (type == "samples") {
    color_scheme = c(brewer.pal(5, "Set1"), brewer.pal(8, "Dark2"), pal_igv("default")(51))
  }
  if (type == "clusters") {
    color_scheme = c( pal_d3("category20")(20), pal_d3("category20b")(20), pal_d3("category20c")(20),pal_igv("default")(51))
  }
  return(color_scheme)
}

### mainText
# savepath = paste0(gsub('/$','',savepath),'/')
# if(!file.exists(savepath)) dir.create(savepath)

DAT = readRDS(dat_rdsPath)
loomDat = readRDS(loomDat_rdsPath)

# change some metadata
loomDat$orig.ident <- gsub("H","H_",substr(colnames(loomDat),start=0,stop=4))
loomDat <- RenameCells(loomDat,new.names= paste0(loomDat$cellname,"_",sapply(strsplit(loomDat$orig.ident,"_"),function(x){x[[2]]})))


# all(paste0(loomDat$orig.ident,"_",loomDat$cellname)==paste0(DAT$orig.ident,"_",DAT$cellname))


loomDat@assays$integrated = DAT@assays$integrated
loomDat@reductions = DAT@reductions
loomDat = AddMetaData(loomDat,DAT@meta.data)
loomDat$param = as.character(loomDat[[param]][,1])

cluster.colors = get_color_scheme('clusters')[1:length(unique(loomDat$param))]
names(cluster.colors) = sort(as.character(unique(loomDat$param)))
cell.colors = sapply(loomDat$param,function(x) {cluster.colors[as.character(x)]})
names(cell.colors) = colnames(loomDat)

loomDat = RunVelocity(object = loomDat, deltaT = 1, kCells = as.numeric(kCells), fit.quantile = as.numeric(fit_quantile),
                              spliced.average = as.numeric(spliced_minAvg),unspliced.average=as.numeric(unspliced_minAvg),
                              reduction = 'umap',group.by='param',ncores=8)
saveRDS(loomDat,file=paste0(savepath,param,"_",'loomDat.rds'))

library(ggplot2)
pdf(paste0(savepath,"cell_velocity_umap_",kCells,'kCells_',n,"_",param,".pdf"),useDingbats=F,width=10,height=10)
DimPlot(loomDat,reduction='umap',group.by='param',cols=cluster.colors) + theme(aspect.ratio=1)
dev.off()

pdf(paste0(savepath,"cell_velocity_",kCells,'kCells_',n,"_",param,".pdf"),useDingbats=F,width=10,height=10)
rs = show.velocity.on.embedding.cor(emb = Embeddings(object = loomDat, reduction = "umap"), 
     vel = Tool(object = loomDat, slot = "RunVelocity"), n = as.numeric(n), scale = "sqrt", 
     cell.colors = ac(x = cell.colors, alpha = 0.5), cex = 0.8, arrow.scale = 3, 
     show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
    do.par = FALSE, cell.border.alpha = 0.1,return.details=TRUE)
dev.off()

save = list(transitionProbability = rs$tp,arrowEstimatesPos = rs$arrows, scale = rs$scale)
saveRDS(save,file=paste0(savepath,param,"_",'velocyto.rds'))








































