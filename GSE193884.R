

# do some data for PS response 
# 2022-4-4
# GSE193884
# download filter data and changed name

library(Seurat)
tmp.dat <- Read10X(data.dir="/public/workspace/lily/PS/Response/GSE193884/")
dat = CreateSeuratObject(counts = tmp.dat)

# do some prepare
tmp_dat = NormalizeData(object = dat)
tmp_dat <- FindVariableFeatures(object = tmp_dat)
# scaling
all.genes <- rownames(x = tmp_dat)
tmp_dat <- ScaleData(object = tmp_dat, features = all.genes)
# PCA
dat <- RunPCA(object = tmp_dat, features = VariableFeatures(object = tmp_dat))


#===================================================================================================================================================================
# read meta data 

metadata <- read.csv("/public/workspace/lily/PS/Response/GSE193884/GSE193884_scRNAseq_multigen_filtered_cells_metadata.csv")
all(colnames(dat)==metadata$cell_name)
dat$sample <- metadata$sample
dat$tumor <- metadata$tumour
dat$gen <- metadata$gen
dat$cell_type <- metadata$cell_type


tumor <- subset(dat,cells=which(dat$cell_type=="malignant"))
saveRDS(tumor,file="/public/workspace/lily/PS/Response/GSE193884/Filter_Tumor_gen.RDS")


#===================================================================================================================================================================
# calculate data 
# ssGSEA

library(Seurat)
dat <- readRDS("/public/workspace/lily/PS/Response/GSE193884/Filter_Tumor_gen.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat[['RNA']]@data),c("p2PSup211216test"),"/public/workspace/lily/PS/",permN=0)
mod <- as.data.frame(mod)
# saveRDS(mod,file="")

dat$PSscore <- mod$p2PSup211216test_norm
dat$group <- paste0(dat$tumor,"_",dat$gen)



# calculate each gen 
library(Seurat)
dat <- readRDS("/public/workspace/lily/PS/Response/GSE193884/Filter_Tumor_gen.RDS")

dat$group <- paste0(dat$tumor,"_",dat$gen)
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
PSscore.res <- list()
samplelist <- as.vector(unique(dat$group))
for(i in 1:length(samplelist)){
	subdat <- subset(dat,cells=which(dat$group==samplelist[i]))
	tmp.mod <- mod.analyze2(as.matrix(subdat[['RNA']]@data),c("p2PSup211216test"),"/public/workspace/lily/PS/",permN=1000)
	tmp.mod <- as.data.frame(tmp.mod)
	tmp.mod$gen <- subdat$gen
	PSscore.res[[i]] <- tmp.mod
	names(PSscore.res)[i] <- samplelist[i]
}


PSscore.res.m <- lapply(PSscore.res,function(x){
	x$flag <- 0
	x$flag[which(x[,3]<0.1)] <- 1 # column 3 is Pvalue and Col 2 is score
	x$gen <- as.vector(x$gen)
	x$flag <- as.numeric(x$flag)
	x
})

# save(PSscore.res,PSscore.res.m,file="/public/workspace/lily/PS/Response/GSE193884/PSscore_Pvalue_0.1.RData")

save(PSscore.res,PSscore.res.m,file="/public/workspace/lily/PS/Response/GSE193884/PSscore_Pvalue_0.1_gen.RData")

res <- lapply(PSscore.res.m,function(x){
	apply(table(x$gen,x$flag),1,function(y){y/sum(y)})
	})



tmp.res <- unlist(res)
tmp.res <- tmp.res[!is.nan(tmp.res)]
tmp.res <- tmp.res[seq(0,length(tmp.res),2)]

gen1 <- tmp.res[c(1,3,9)] # delete 142
gen2 <- tmp.res[c(2,4,10)]  # delete 142
gen3 <- tmp.res[c(5,8,11)]


# plot result 
# gen1 <- sapply(res,function(x){x[2,1]})
# gen2 <- sapply(res,function(x){x[2,2]})
# gen3 <- sapply(res,function(x){x[2,3]})[is.finite(sapply(res,function(x){x[2,3]}))] # have one NaN

pdata <- data.frame(percent = c(gen1,gen2,gen3),group = c(rep("gen1",length(gen1)),rep("gen2",length(gen2)),rep("gen3",length(gen3))))
# pdf("/public/workspace/lily/PS/Response/GSE193884/PSscore_Pvalue_0.1_gen.pdf",useDingbats=F)
boxplot(percent~group,data=pdata,names=c("gen1","gen2","gen3"),main="Pvalue <0.1")
beeswarm::beeswarm(percent~group,data=pdata,pch = 16, add=T)



dev.off()





















