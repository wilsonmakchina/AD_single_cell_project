# roc curve
path <- '/share/home/duanshumin/wilsonmak/analysis/AD_microglia/GSE103334/'

pkgs <- list("dplyr","DMwR","caret", "glmnet", "doParallel", "foreach", "pROC", "stabs","parallel","tidyr")

source("/share/home/duanshumin/wilsonmak/analysis/AD_microglia/machine_learning_functions.R")

lapply(pkgs, require, character.only = T)

load('/share/home/duanshumin/wilsonmak/analysis/AD_microglia/GSE103334/machine_learning/20181104/GSE103334_stabs_selected_genes.Rdata')

set.seed(2018)
randomly_select_genes1 <- sample(colnames(learnCV), 47)
randomly_select_genes1 <- which(colnames(learnCV[-1]) %in% randomly_select_genes1)

set.seed(2019)
randomly_select_genes2 <- sample(colnames(learnCV), 47)
randomly_select_genes2 <- which(colnames(learnCV[-1]) %in% randomly_select_genes2)

set.seed(2020)
randomly_select_genes3 <- sample(colnames(learnCV), 47)
randomly_select_genes3 <- which(colnames(learnCV[-1]) %in% randomly_select_genes3)



acc <- NULL
ROC <- NULL
p <- NULL
auc <- NULL
selected_genes <- NULL
selected_genes[[1]] <- c(1:ncol(testCV[-1]))
selected_genes[[2]] <- stab$selected
selected_genes[[3]] <- randomly_select_genes1
selected_genes[[4]] <- randomly_select_genes2
selected_genes[[5]] <- randomly_select_genes3


control <- trainControl(method = 'none', classProbs = TRUE, summaryFunction = twoClassSummary)

for (i in 1:5) {
	svm <- train(learnCV[-1][selected_genes[[i]]], learnCV[,1], method = "svmLinear",trControl = control, metric='ROC')
	predictY <- predict(svm, testCV[-1][selected_genes[[i]]], type = "raw")
	acc[i] <- mean(testCV$y == predictY)
	svm.probs <- predict(svm, testCV[-1][selected_genes[[i]]], type = "prob")
	ROC[[i]] <- roc(response = testCV$y, predictor = svm.probs$cluster_5, levels = levels(testCV$y))
	auc[i] <- ROC[[i]]$auc
	if (i <= 2) {
		p[[i]] <- ggroc(list(roc_curve=ROC[[i]]),legacy.axes = TRUE) 
	}
}

p[[3]] <- ggroc(list(random_selected_1=ROC[[3]], random_selected_2=ROC[[4]], random_selected_3=ROC[[5]]),legacy.axes = TRUE)


p[[4]] <- ggroc(list(all=ROC[[1]], ML_selected=ROC[[2]], random_selected_1=ROC[[3]], random_selected_2=ROC[[4]], random_selected_3=ROC[[5]]),legacy.axes = TRUE)


# pdf('roc_curve.pdf')
# print(p)
# dev.off()
save(p, acc, auc, file='performance_GSE103334.Rdata')
