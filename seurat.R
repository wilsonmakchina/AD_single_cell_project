library(Seurat)
library(dplyr)
library(tidyr)
library(Matrix)

path = '/share/home/duanshumin/wilsonmak/analysis/AD_microglia/GSE98969/'



cell_list <- read.table(paste0(path,'MG_list.txt'), sep = '\t', stringsAsFactors = F, head = 1)
# cell_list <- cell_list[cell_list$seuratCluster %in% c('2: Resting MG', '4: DAM'),] 


MG <- read.csv(paste0(path,'all_cell_raw_counts.csv'),stringsAsFactors = F, row.names=1, head=1)

MG <- MG[,cell_list$cellID]

MG <- CreateSeuratObject(raw.data = MG)

rownames(cell_list) = cell_list$cellID
cell_type <- cell_list[4]
colnames(cell_type) <- 'cell_type'
a <- cell_type %>% separate(cell_type, c("cluster", "type"), ":")
a$cluster <- as.numeric(a$cluster)
a <- a[order(a$cluster),]
cell_type1 <- merge(a, cell_type, by = 0)
rownames(cell_type1) <- cell_type1$Row.names
cell_type1 <- cell_type1[order(cell_type1$cluster),]
cell_type1$cell_type <- factor(cell_type1$cell_type, levels=unique(cell_type1$cell_type))
MG <- AddMetaData(object = MG, metadata = cell_type1, col.name = "cell_type")



MG <- NormalizeData(object = MG, normalization.method = "LogNormalize",
                      scale.factor = 10000)

MG <- FindVariableGenes(object = MG, mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = MG@var.genes)
MG <- ScaleData(object = MG)
MG <- RunPCA(object = MG, pc.genes = MG@var.genes, do.print = TRUE, pcs.print = 1:5, 
    genes.print = 5)


PCElbowPlot(object = MG)


MG <- RunTSNE(object = MG, dims.use = 1:15, do.fast = T, perplexity=30)




p <- TSNEPlot(object = MG, group.by = "strain", do.return = TRUE, pt.size = 1)



pdf(paste0('tsne_cell_type_for_legend.pdf'), width = 5, height = 5)
TSNEPlot(object = MG, group.by = "cell_type", do.return = TRUE, pt.size = 0.1)
dev.off()





MG@ident <- as.factor(cell_type1$type)
names(MG@ident) <- rownames(cell_type1)


gene_name = '4632428N05Rik'
pdf(paste0(gene_name, '_vlnPlot.pdf'))
VlnPlot(object = MG, features.plot = c(gene_name), ident.include = c(' Resting MG', ' DAM') )
dev.off()

saveRDS(MG, file = "MG.rds")

