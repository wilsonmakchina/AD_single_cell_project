library(Seurat)
library(dplyr)
library(Matrix)

## figure margin
par(mar=rep(2,4))


path = '/share/home/duanshumin/wilsonmak/analysis/AD_microglia/GSE103334/'


MG_data <- read.table(paste0(path,'GSE103334_FPKM_CKP25_TOPHAT.txt'), header = T, stringsAsFactors = F)
MG_data <- tibble::column_to_rownames(MG_data, var= 'GENES')



cell_list <- read.table(paste0(path,'1685_Microglia.txt'), stringsAsFactors = F)
MG <- MG_data[,c(cell_list$V1)]
gene_filter_setting <- rowSums(MG > 0) >= 30
MG_filter_genes <- row.names(MG[gene_filter_setting,])
length(MG_filter_genes)
MG <- CreateSeuratObject(raw.data = MG)

metadata <- read.table(paste0(path,'metadata.txt'), head=F, row.name=4, stringsAsFactors=FALSE)
colnames(metadata) <- c('strain','time_point','biological_replicates')
MG <- AddMetaData(object = MG, metadata = metadata[1], col.name = "strain")

strain_time_point <- as.data.frame(factor(unlist(
    lapply(colnames(MG@data),
           function(x) strsplit(as.character(x), '_m')[[1]][1]))), row.names = colnames(MG@data))
colnames(strain_time_point) <- "strain_time_point"
MG <- AddMetaData(object = MG, metadata = strain_time_point, col.name = "strain_time_point")

MG <- AddMetaData(object = MG, metadata = metadata[3], col.name = "biological_replicates")



MG <- NormalizeData(object = MG, normalization.method = "LogNormalize",
                      scale.factor = 10000)

a <- MG@data
rownames(a)
a <- a[rownames(MG@scale.data),]
a <- as.matrix(a)
rownames(a)
write.csv(a, file = 'normalized_data.csv')



MG@var.genes = MG_filter_genes
length(x = MG@var.genes)

MG <- ScaleData(object = MG, genes.use = MG_filter_genes)#, vars.to.regress = "biological_replicates")

## standard deviation of PCs
MG <- RunPCA(object = MG, pc.genes = MG@var.genes, do.print = T, pcs.print = 1:5,
               genes.print = 5)
pdf('PCA.pdf', width = 5, height = 4)
PCElbowPlot(object = MG, num.pc = 20)
dev.off()


MG <- RunTSNE(object = MG, dims.use = 1:10, do.fast = F, perplexity=30)
MG <- DBClustDimension(MG, reduction.use = "tsne", G.use = 5)

pdf(paste0('tsne.pdf'), width = 5, height = 3)
TSNEPlot(object = MG, pt.size = 1)
dev.off()

## CK and CKp-25 stne
pdf(paste0('tsne_strain.pdf'), width = 5, height = 3)
TSNEPlot(object = MG, group.by = "strain", do.return = TRUE, pt.size = 1)
dev.off()
## time point: CK and CKp-25 stne
p <- NULL
for(i in 1:8){
    strain_time_point <- levels(factor(unlist(
        lapply(colnames(MG@data),
               function(x) strsplit(as.character(x), '_m')[[1]][1]))))
    # pdf(paste0('tsne_',strain_time_point[i], '.pdf'), width = 5, height = 3)
    colors <- c(rep('gray',8))
    colors[i] <- 2
    p[[i]] <- TSNEPlot(object = MG, group.by = "strain_time_point", 
             colors.use = colors, do.return = TRUE, pt.size = 1)
    # print(p)
    # dev.off()
}
library(cowplot)
plot = plot_grid(plotlist = p,align = 'v', ncol = 2)
save_plot("t-sne_time_point.pdf", plot, base_height = 11,base_width = 10,
          base_aspect_ratio = 2 # make room for figure legend
)

## biological replicates stne
TSNEPlot(object = MG, group.by = "biological_replicates", 
         do.return = TRUE, pt.size = 2)

a <- MG@cell.names[(MG@meta.data$DBclust.ident == 3 | MG@meta.data$DBclust.ident == 6)]
a <- MG@cell.names[MG@meta.data$DBclust.ident == 9]

write.table(a, 'cluster9_need_to_be_remove.txt', row.names = F, col.names = F, quote = F)
length(a)
length(grep('_6w_', a))
a <- grep('CK_', colnames(MG@raw.data))
length(a)
b <- grep('CKp25_',colnames(MG@data))
length(b)
save(MG, file = "~/Desktop/project/inbox/GSE103334_single_cell_ML_analysis/2018.1.16/MG_filter_at_least_30_cell_contain_each_gene_11_clusters.Rda")


scale_data <- as.data.frame(MG@scale.data)
scale_data <- scale_data[,(F == colnames(scale_data) %in% a)]
ncol(scale_data)
class(scale_data)
write.csv(scale_data, 'scale_data_for_supervised_learning.csv')


a <- MG@cell.names[(MG@meta.data$DBclust.ident == 9)]
library(pheatmap)
optimal_svm_genes_6w <- as.list(read.csv('optimal_ranked_genes_90.csv')[3])
rownames(MG_data)
data_heatmap <- MG@raw.data[rownames(MG@raw.data) %in% optimal_svm_genes_6w$X1,]
data_heatmap <- data_heatmap[,(F == colnames(data_heatmap) %in% a)]
list_6w <- grep('_6w', colnames(data_heatmap))
data_heatmap <- data_heatmap[,list_6w]


pdf('heatmap_90_genes.pdf', width = 15, height = 17)
pheatmap(log2(data_heatmap+1), scale = "row", border_color= 'NA', fontsize = 15, cluster_cols = F,
         color=colorRampPalette(c("#008800","#008800","green","green",  "black","yellow", "yellow","#DDAA00","#DDAA00"))(500))
dev.off()
library(ggplot2)
