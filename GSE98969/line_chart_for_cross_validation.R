library(ggplot2)
library(foreach)

setwd("/Volumes/data_wilson/OneDrive/OneDrive - zju.edu.cn/GSE98969/20181103")

data_name = 'GSE98969'
a <- system('ls CV*.csv', intern = T)

stats <- foreach(file_name = a, .combine = rbind) %do% {
    df <- read.csv(file=file_name)
    file_name = rev(strsplit(file_name,'_')[[1]])[1]
    if (file_name == 'lasso.csv') {file_name='1.csv'}
    file_name <- gsub('.csv','',file_name)
    df[4] <- file_name
    df
}
# write.csv(stats, 'CV_result_lasso_or_net.csv', row.names=F)



gg <- ggline(stats, x = "cutoff", y = "Average", color = 'V4',size=1, palette = "Set1") +
    labs(title=data_name,x=expression(pi[thr]), y="Accuracy", color = expression(paste(alpha,' value'))) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.margin = unit(c(1,5,1,5), "mm")) +
    theme(text = element_text(size = 15)) 
    
    
    

setwd("/Volumes/data_wilson/OneDrive/OneDrive - zju.edu.cn/GSE103334_single_cell_ML_analysis/20181103")

data_name = 'GSE103334'
a <- system('ls CV*.csv', intern = T)

stats <- foreach(file_name = a, .combine = rbind) %do% {
    df <- read.csv(file=file_name)
    file_name = rev(strsplit(file_name,'_')[[1]])[1]
    if (file_name == 'lasso.csv') {file_name='1.csv'}
    file_name <- gsub('.csv','',file_name)
    df[4] <- file_name
    df
}


gg2 <- ggline(stats, x = "cutoff", y = "Average", color = 'V4',size=1, palette = "Set1") +
    labs(title=data_name,x=expression(pi[thr]), y="Accuracy", color = expression(paste(alpha,' value'))) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.margin = unit(c(1,5,1,5), "mm")) +
    theme(text = element_text(size = 15))

p <- ggarrange(gg, gg2, labels = c("A", "B"),  ncol = 2, nrow = 1)

pdf('cross_validation.pdf', width = 10, height = 3)
print(p)
dev.off()
