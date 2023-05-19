
# R 
# 202.195.187.5
library(estimate)
# calculate Brain metastasis sample purity 
#===============================================================================================================================================
# dat <- read.table("/public/workspace/wulx/tst/CGGA.mRNAseq_325.RSEM-genes.20200506.txt",header=T)
# write.table(dat,file="GSE14108_exp.txt",sep="\t",row.names=T,col.names=T,quote=F)
filterCommonGenes(input.f="/public/workspace/wulx/tst/CGGA.mRNAseq_325.RSEM-genes.20200506.txt",output.f="~/tmp/CGGA_325.gct",id="GeneSymbol")
estimateScore("~/tmp/CGGA_325.gct", "~/tmp/CGGA_325_estimate_score.gct", platform="affymetrix")



filterCommonGenes(input.f="/public/workspace/wulx/tst/CGGA.mRNA_array_301_gene_level.20200506.txt",output.f="~/tmp/CGGA_301.gct",id="GeneSymbol")
estimateScore("~/tmp/CGGA_301.gct", "~/tmp/CGGA_301_estimate_score.gct", platform="affymetrix")


















