library(Seurat)
library(dplyr)
library(Matrix)

path = '/share/home/duanshumin/wilsonmak/analysis/AD_microglia/GSE98969/'

# a <- list.files()
# dat <- read.table(a[1], row.names = 1, head = 1)
# for (i in 2:97) {
# 	df <- read.table(a[i], row.names = 1, head = 1)
# 	dat <- cbind(dat, df)
# }


cell_list <- read.table(paste0(path,'MG_list.txt'), sep = '\t', stringsAsFactors = F, head = 1)
cell_list <- cell_list[cell_list$seuratCluster %in% c('2: Resting MG', '4: DAM'),] 

# MG <- dat[,cell_list$cellID]

MG <- read.csv('resting_MG_and_DAM.csv',stringsAsFactors = F, row.names=1, head=1)

gene_filter_setting <- rowSums(MG > 0) >= 30
MG_filter_genes <- row.names(MG[gene_filter_setting,])
length(MG_filter_genes)
MG <- CreateSeuratObject(raw.data = MG)
cell_type <- cell_list[4]
colnames(cell_type) <- 'cell_type'
MG <- AddMetaData(object = MG, metadata = cell_type, col.name = "cell_type")

MG <- NormalizeData(object = MG, normalization.method = "LogNormalize",
                      scale.factor = 10000)

MG <- RunTSNE(object = MG, dims.use = 1:5, do.fast = F, perplexity=30)
MG <- DBClustDimension(MG, reduction.use = "tsne", G.use = 5)



p <- TSNEPlot(object = MG, group.by = "strain", do.return = TRUE, pt.size = 1)


pdf(paste0('tsne_cell_type.pdf'), width = 5, height = 3)
TSNEPlot(object = MG, group.by = "cell_type", do.return = TRUE, pt.size = 0.1)
dev.off()
