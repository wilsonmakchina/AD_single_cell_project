library(Seurat)
library(dplyr)
library(tidyr)
library(Matrix)
library(RColorBrewer)


path = '/share/home/duanshumin/wilsonmak/analysis/AD_microglia/GSE98969/'



cell_list <- read.table(paste0(path,'MG_list.txt'), sep = '\t', stringsAsFactors = F, head = 1)
# cell_list <- cell_list[cell_list$seuratCluster %in% c('2: Resting MG', '4: DAM'),] 


MG <- read.csv(paste0(path,'all_cell_raw_counts.csv'),stringsAsFactors = F, row.names=1, head=1)

MG <- MG[,cell_list$cellID]



rownames(cell_list) = cell_list$cellID
cell_type <- cell_list[4]
colnames(cell_type) <- 'cell_type'
a <- cell_type %>% separate(cell_type, c("cluster", "type"), ": ")
a$cluster <- as.numeric(a$cluster)
a <- a[order(a$cluster),]
cell_type1 <- merge(a, cell_type, by = 0)
rownames(cell_type1) <- cell_type1$Row.names
cell_type1 <- cell_type1[order(cell_type1$cluster),]
cell_type1$cell_type <- factor(cell_type1$cell_type, levels=unique(cell_type1$cell_type))


MG <- CreateSeuratObject(raw.data = MG)
mito.genes <- grep(pattern = "^mt-", x = rownames(x = MG@raw.data), value = TRUE)
percent.mito <- Matrix::colSums(MG@raw.data[mito.genes,])/Matrix::colSums(MG@raw.data)

# Calculate ERCC abundances on the raw counts before creating a Seurat object
ERCC.index <- grep(pattern = "^ERCC-", x = rownames(MG@raw.data), value = FALSE) # Select row indices and not ERCC names 
percent.ERCC <- Matrix::colSums(MG@raw.data[ERCC.index,])/Matrix::colSums(MG@raw.data)

# Remove ERCC from MG

MG <- MG@raw.data[-ERCC.index,]

# Create Seurat object, and add percent.ERCC to object@meta.data in the percent.ERCC column
MG <- CreateSeuratObject(raw.data = MG)



MG <- AddMetaData(object = MG, metadata = cell_type1, col.name = "cell_type")
MG <- AddMetaData(object = MG, metadata = percent.mito, col.name = "percent.mito")
MG <- AddMetaData(object = MG, metadata = percent.ERCC, col.name = "percent.ERCC")


MG <- NormalizeData(object = MG, normalization.method = "LogNormalize",
                      scale.factor = 10000)


MG <- ScaleData(object = MG, do.scale = T, do.center = T, vars.to.regress = c("nUMI", "percent.mito", "percent.ERCC"))

# MG <- ScaleData(object = MG)
MG <- FindVariableGenes(object = MG, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = MG@var.genes)
MG <- RunPCA(object = MG, pc.genes = MG@var.genes, do.print = TRUE, pcs.print = 1:5, 
    genes.print = 5)


# PCElbowPlot(object = MG)


MG <- RunTSNE(object = MG, dims.use = 1:15, do.fast = T, perplexity=30)




p <- TSNEPlot(object = MG, group.by = "strain", do.return = TRUE, pt.size = 1)

strain <- read.csv(paste0(path,'strain_information.csv'),stringsAsFactors = F, head = 1, row.names=1)

MG <- AddMetaData(object = MG, metadata = strain, col.name = "strain")


pdf(paste0('tsne_cell_type.pdf'), width = 5, height = 3)
TSNEPlot(object = MG,  group.by='Mouse_ID',do.return = TRUE, pt.size = 0.1)
dev.off()



MG@ident <- as.factor(cell_type1$type)
names(MG@ident) <- rownames(cell_type1)

color = colorRampPalette(c('grey', brewer.pal(9,"Blues"), rev(brewer.pal(9,'YlOrBr')), brewer.pal(9,"YlOrRd")))(100)
# color = colorRampPalette(brewer.pal(12,'Paired'))(100)
list = names(MG@ident)[MG@ident == ' DAM' | MG@ident == ' Resting MG']
gene_name = 'A130022J15Rik'
pdf(paste0(gene_name, '_featurePlot.pdf'))
p <- FeaturePlot(object = MG, features.plot = c(gene_name), reduction.use = "tsne",cells.use = list, cols.use = color)
print(p)
dev.off()

for (i in colnames(mydata)[stab$selected]) {
gene_name = i
pdf(paste0(gene_name, '_vlnPlot.pdf'))
p <- VlnPlot(object = MG, features.plot = c(gene_name), ident.include = c('Resting MG', 'DAM') )
print(p)
dev.off()
}
saveRDS(MG, file = "MG_regressed_out.rds")

