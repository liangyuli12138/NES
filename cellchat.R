
# 2021-5-29
# Rscript 
# run analysis for Tumor and Macrophage interaction
###################################################################################################
library(Seurat)
library(CellChat)

dat1 <- readRDS("/public/workspace/wulx/missions/Macrophage2tumor/rd1mm.RDS")
# Create a cellchat obj 
data.input <- GetAssayData(dat1, assay = "RNA", slot = "data") # normalized data matrix
meta <- data.frame(dat1@meta.data) # create a dataframe of the cell labels

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
# change ident and add meta information
# cellchat <- addMeta(cellchat, meta = meta, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "celltype") # set "labels" as default cell identity
levels(cellchat@idents)


# database
###################################################################################################
CellChatDB <- CellChatDB.human
# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB) 
# use all DB 
CellChatDB.use <- CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use

# 
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel


#####################################################################################################
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# project gene expression data onto PPI network (should use)
cellchat <- projectData(cellchat, PPI.human)

# Calculate cell cell communication 
####################################################################################################
cellchat <- computeCommunProb(cellchat) # 1W cells 5min
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)




cellchat <- computeCommunProbPathway(cellchat)


# Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
#par(mfrow = c(1,2), xpd=TRUE)



# interaction pair 
#netVisual_bubble(cellchat, remove.isolate = T)

# interaction pair for bubble plot 
tmp.dat <- subsetCommunication(cellchat)
tmp.dat$source.target <- paste0(tmp.dat$source," -> ",tmp.dat$target)
df <- tmp.dat[which(tmp.dat$prob>0.01&tmp.dat$pval<0.01),]



pdf("/public/workspace/lily/wulx/dat1.pdf",useDingbats=F,height=10,width=10)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

ggplot(df, aes(x = source.target, y = interaction_name_2,
	color = prob, size = pval)) + geom_point(pch = 16) +
	theme_linedraw() + theme(panel.grid.major = element_blank()) +
	theme(axis.text.x = element_text(angle = 90, hjust = 1,
	vjust = 0.5), axis.title.x = element_blank(),
	axis.title.y = element_blank()) + scale_x_discrete(position = "bottom")


dev.off()


saveRDS(cellchat,file="/public/workspace/lily/wulx/dat1_cellchat.RDS")



######################################################################################################################################
# dat2 

library(Seurat)
library(CellChat)

dat2 <- readRDS("/public/workspace/wulx/missions/Macrophage2tumor/rd2mm.RDS")
# Create a cellchat obj 
data.input <- GetAssayData(dat2, assay = "RNA", slot = "data") # normalized data matrix
meta <- data.frame(dat2@meta.data) # create a dataframe of the cell labels

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
# change ident and add meta information
# cellchat <- addMeta(cellchat, meta = meta, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "celltype") # set "labels" as default cell identity
levels(cellchat@idents)


# database
###################################################################################################
CellChatDB <- CellChatDB.human
# use all DB 
CellChatDB.use <- CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel
#####################################################################################################
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# project gene expression data onto PPI network (should use)
cellchat <- projectData(cellchat, PPI.human)

# Calculate cell cell communication 
####################################################################################################
cellchat <- computeCommunProb(cellchat) # 1W cells 5min
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)




cellchat <- computeCommunProbPathway(cellchat)


# Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
#par(mfrow = c(1,2), xpd=TRUE)



# interaction pair 
#netVisual_bubble(cellchat, remove.isolate = T)

# interaction pair for bubble plot 
tmp.dat <- subsetCommunication(cellchat)
tmp.dat$source.target <- paste0(tmp.dat$source," -> ",tmp.dat$target)
df <- tmp.dat[which(tmp.dat$prob>0.01&tmp.dat$pval<0.01),]



pdf("/public/workspace/lily/wulx/dat2.pdf",useDingbats=F,height=10,width=10)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

ggplot(df, aes(x = source.target, y = interaction_name_2,
	color = prob, size = pval)) + geom_point(pch = 16) +
	theme_linedraw() + theme(panel.grid.major = element_blank()) +
	theme(axis.text.x = element_text(angle = 90, hjust = 1,
	vjust = 0.5), axis.title.x = element_blank(),
	axis.title.y = element_blank()) + scale_x_discrete(position = "bottom")

dev.off()


saveRDS(cellchat,file="/public/workspace/lily/wulx/dat2_cellchat.RDS")


















































