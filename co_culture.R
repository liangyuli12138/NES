
# 2021-9-21
# analysis for co-culture 
#==================================================================================================================================================
tmp <- read.table("/public/workspace/lily/wulx/coculture2_3v3.txt")

dat <- as.matrix(tmp[,c(2:7)])

source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
modlist <- gsub("\\.mod","",dir("/public/workspace/lily/MOD_file/HALLMARK/"))
mod <- mod.analyze2(as.matrix(dat),modlist,"/public/workspace/lily/MOD_file/HALLMARK/",permN=0)
mod <- as.data.frame(mod)


mod.f <- t(mod[,c(51:100)])
# get some pathway 
dat.score <- mod.f[c(29,44,22,5),]


library(pheatmap)
# add a annotation 
ann <- data.frame(row.names=colnames(dat.score),group=c(rep("monoculture",3),rep("co_culture",3)),
	time=c("24h","48h","72h","24h","48h","72h"))
annoCol<-list(group=c(monoculture="#393B94", co_culture="#F8C26B"))

pdf("/public/workspace/lily/wulx/culture_res.pdf",useDingbats=F)
pheatmap(dat.score,cluster_cols=F,annotation_col=ann,gaps_col=3,cellwidth=20,cellheight=20,
	color = colorRampPalette(colors = c("steelblue","white","red"))(100),annotation_colors= annoCol)
dev.off()





# extravasion signature plot
tmp <- read.table("/public/workspace/lily/wulx/coculture2_3v3.txt")

dat <- as.matrix(tmp[,c(2:7)])

source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
meta_gene <- c("SNAI1","SNAI2","TWIST1","ZEB1","ZEB2",'TGFB1','CSF1','EGF','VEGFA','PTGS2','EREG','ANGPT2','MMP1','MMP2','MMP3','MMP10')
mod.generate(meta_gene,"meta_gene",out="/public/workspace/lily/MOD_file/meta_gene.mod") # make a mod file 

#modlist <- gsub("\\.mod","",dir("/public/workspace/lily/MOD_file/HALLMARK/"))
mod <- mod.analyze2(as.matrix(dat),"meta_gene","/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)

dat.f <- dat[meta_gene,]
 all(rownames(mod)==colnames(dat.f))
library(pheatmap)
ann <- data.frame(row.names=colnames(dat.f),group=c(rep("monoculture",3),rep("co_culture",3)),
	time=c("24h","48h","72h","24h","48h","72h"),metastasis_sig=mod[,2]) # this NES is extravasion score # 2021-9-21

pdf("/public/workspace/lily/wulx/culture_extravasion_res.pdf",useDingbats=F)
pheatmap(dat.f,scale="row",cluster_cols=F,cluster_rows=F,annotation_col=ann,gaps_col=3,cellwidth=20,cellheight=20,
	color = colorRampPalette(colors = c("steelblue","white","red"))(100))
dev.off()








#=====================================================================================================================================================
# data for PS 
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
GSE27523 = read.table('/public/workspace/lily/wulx/GSE27523_expr.txt',header=T,sep='\t')
gds = GEOquery::getGEO(filename='/public/workspace/lily/wulx/GSE27523_family.soft')
map <- gds@gpls[[1]]@dataTable@table
map$Gene = sapply(as.vector(map$gene_assignment),function(x){
  unlist(strsplit(unlist(strsplit(x,' // ')[[1]][2]),' /// ')[[1]][1])
})
ref = map[,c('ID','Gene')]
ref <- ref[which(!is.na(ref$Gene)),]
colnames(ref) <- c('ID_REF','Gene')
GSE27523$ID_REF = rownames(GSE27523)
GSE27523.expr <- merge(by='ID_REF',GSE27523,ref)[,-1]
GSE27523.expr <- aggregate(.~Gene,data=GSE27523.expr,mean)
rownames(GSE27523.expr) <- GSE27523.expr$Gene
GSE27523.expr <- GSE27523.expr[,-1]
rs = as.data.frame(t(mod.analyze(GSE27523.expr,c('p2PSup0714',"HALLMARK_HYPOXIA"),mod.dir="/public/workspace/lily/PS/Final_716/")))
info = read.table('/public/workspace/lily/wulx/GSE27523_pheno.txt',header=T,row.names = 1,sep='\t')
info$type=0
info[1:3,]$type= 1 # control
info[4:6,]$type= 2 # H1A
info[7:9,]$type= 3 # H2A


info$ps = rs$p2PSup0714
info$hyp = rs$HALLMARK_HYPOXIA
pdf("/public/workspace/lily/wulx/HIFA_PS.pdf",useDingbats=F,width=5)
boxplot(ps~type,info[-c(10:12),],boxwex=0.5,names=c("control","H1A","H2A"),las=2)
beeswarm::beeswarm(ps~type,info[-c(10:12),],cex=2,add=T)
dev.off()



wilcox.test(ps~type,info[c(1:6),])
wilcox.test(ps~type,info[c(1:3,7:9),])



#========================================================================================
info$FOSL2 <- as.numeric(as.vector(GSE27523.expr["FOSL2",]))

pdf("/public/workspace/lily/wulx/HIFA_FOSL2.pdf",useDingbats=F,width=5)
boxplot(FOSL2~type,info[-c(10:12),],boxwex=0.5,names=c("control","H1A","H2A"),las=2)
beeswarm::beeswarm(FOSL2~type,info[-c(10:12),],cex=2,add=T)
dev.off()


wilcox.test(FOSL2~type,info[c(1:6),])
wilcox.test(FOSL2~type,info[c(1:3,7:9),])














#=================================================================================================================================
# lesion1 calculate 
library(Seurat)
dat <- readRDS("/public/workspace/lily/PS/data/hbrd001.rds")
data <- dat[["RNA"]]@data

dat$CCL2 <- ifelse(data["CCL2",]>0,"Y","N")
dat$CCR10 <- ifelse(data["CCR10",]>0,"Y","N")

table(dat$hbmarker,dat$CCL2)
































