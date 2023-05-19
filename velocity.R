
# 2021-10-8
# this program is used for analysis velocity 
#==============================================================================================================================================


# 1. lesion1 and lesion2 sample 
# run in 202.195.187.3
# use velocyte to get loom files
conda activate velocity
module load  samtools-1.9

velocyto run10x -m /public/workspace/lily/REF/hg19_rmsk.gtf \
    -@ 10 --samtools-memory 2000 \
    ~/PS/Velocyto/RD0012/RD-20180817-001-SR18271/ \
    ~/REF/refdata-cellranger-hg19-1.2.0/genes/genes.gtf


conda activate velocity
module load  samtools-1.9

velocyto run10x -m /public/workspace/lily/REF/hg19_rmsk.gtf \
    -@ 10 --samtools-memory 2000 \
    ~/PS/Velocyto/RD0012/RD-20180817-002-SR18271/ \
    ~/REF/refdata-cellranger-hg19-1.2.0/genes/genes.gtf



# Seurat object use to tumor loom 
module load hdf5-1.8.13
library(Seurat)
library(loomR)
dat1 <- readRDS("/public/workspace/lily/PS/data/hbrd001.rds")
tumor1 <- subset(dat1,cells=which(dat1$hbmarker=="Tumor cell"))

as.loom(tumor1, assay="RNA", filename = "/public/workspace/lily/PS/Multifocal/Tumor_loom/lesion1_tumor.loom")


dat2 <- readRDS("/public/workspace/lily/PS/data/hbrd002.rds")
tumor2 <- subset(dat2,cells=which(dat2$hbmarker=="Tumor cell"))

as.loom(tumor2, assay="RNA", filename = "/public/workspace/lily/PS/Multifocal/Tumor_loom/lesion2_tumor.loom")











# 2. P912_left and P912_right 
#==================================================================================================================================================
library(Seurat)
library(loomR)
P912L <- readRDS("/public/workspace/lily/GBM/all_cell_celltype/P912_left0514.RDS")
tumor.l <- subset(P912L,cells=which(P912L$type=="malignant"))
as.loom(tumor.l, assay="RNA", filename = "/public/workspace/lily/PS/Multifocal/Tumor_loom/P912L_tumor.loom")


P912R <- readRDS("/public/workspace/lily/GBM/all_cell_celltype/P912_right0514.RDS")
tumor.R <- subset(P912R,cells=which(P912R$type=="malignant"))
as.loom(tumor.R, assay="RNA", filename = "/public/workspace/lily/PS/Multifocal/Tumor_loom/P912R_tumor.loom")



conda activate velocity
module load  samtools-1.9

velocyto run10x -m /public/workspace/lily/REF/hg19_rmsk.gtf \
    -@ 10 \
    ~/PS/Velocyto/P912_left0514/ \
    ~/REF/refdata-cellranger-hg19-1.2.0/genes/genes.gtf



conda activate velocity
module load  samtools-1.9

velocyto run10x -m /public/workspace/lily/REF/hg19_rmsk.gtf \
    -@ 10 \
    ~/PS/Velocyto/P912_right0514/ \
    ~/REF/refdata-cellranger-hg19-1.2.0/genes/genes.gtf









# 3. P673 and P689 sample 
#===================================================================================================================================================
library(Seurat)
library(loomR)
P673 <- readRDS("/public/workspace/lily/GBM/all_cell_celltype/P673_jc.RDS")
tumor.P673 <- subset(P673,cells=which(P673$type=="malignant"))
as.loom(tumor.P673, assay="RNA", filename = "/public/workspace/lily/PS/Multifocal/Tumor_loom/P673_tumor.loom")


P689 <- readRDS("/public/workspace/lily/GBM/all_cell_celltype/P689_2.RDS")
tumor.P689 <- subset(P689,cells=which(P689$type=="malignant"))
as.loom(tumor.P689, assay="RNA", filename = "/public/workspace/lily/PS/Multifocal/Tumor_loom/P689_tumor.loom")




# this 
conda activate velocyto
module load  samtools-1.9

velocyto run10x -m /public/workspace/lily/REF/hg19_rmsk.gtf \
    -@ 30 --samtools-memory 2000 \
    /public/workspace/lily/PS/Multifocal/Result/P673_jc/ \
    ~/REF/INDEX-hg19/refdata-cellranger-hg19-1.2.0/genes/genes.gtf



conda activate velocyto
module load  samtools-1.9

velocyto run10x -m /public/workspace/lily/REF/hg19_rmsk.gtf \
    -@ 30 --samtools-memory 2000 \
    /public/workspace/lily/PS/Multifocal/Result/P689_2/P689_2/ \
    ~/REF/INDEX-hg19/refdata-cellranger-hg19-1.2.0/genes/genes.gtf









# 4. New paired LCBM
#===============================================================================================================================================
# this was run in 202.195.187.3
conda activate velocity
module load  samtools-1.9

velocyto run10x -m /public/workspace/lily/REF/hg19_rmsk.gtf \
    -@ 20 --samtools-memory 2000 \
    /public/workspace/lily/PS/10XT1H3/ \
     ~/REF/refdata-cellranger-hg19-1.2.0/genes/genes.gtf



conda activate velocity
module load  samtools-1.9

velocyto run10x -m /public/workspace/lily/REF/hg19_rmsk.gtf \
    -@ 20 --samtools-memory 2000 \
    /public/workspace/lily/PS/10XT1H4/ \
     ~/REF/refdata-cellranger-hg19-1.2.0/genes/genes.gtf






















