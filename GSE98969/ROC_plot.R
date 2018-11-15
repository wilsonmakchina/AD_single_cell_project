library(ggplot2)
library(ggpubr)

load('/Volumes/data_wilson/OneDrive/OneDrive - zju.edu.cn/GSE98969/20181115/performance_GSE98969.Rdata')

for (i in 1:4) {
    p[[i]] <- p[[i]] + theme_bw() +
        scale_color_brewer(palette = "Set1") +
        labs(color='') +
        theme(legend.position = c(0.7, 0.4)) +
        theme(panel.border = element_rect(colour = "black"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
        geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color='gray', linetype='dashed')
}

gg <- ggarrange(p[[1]], p[[2]], p[[3]], p[[4]], labels = c("C", "D", "E", "F"),  ncol = 2, nrow = 2)
gg

pdf('roc.pdf', width = 5, height = 5)
print(gg)
dev.off()


acc[1]
auc[1]

acc[2]
auc[2]


random_acc <- c(acc[3],acc[4],acc[5])
mean(random_acc)
sd(random_acc)

random_auc <- c(auc[3],auc[4],auc[5])
mean(random_auc)
sd(random_auc)
