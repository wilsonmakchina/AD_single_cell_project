path <- '/share/home/duanshumin/wilsonmak/analysis/AD_microglia/GSE98969/'

pkgs <- list("dplyr","DMwR","caret", "glmnet", "doParallel", "foreach", "pROC", "stabs","parallel","tidyr")
# source(paste0(path,'machine_learning_functions.R'))
source('~/wilsonmak/analysis/AD_microglia/GSE103334/machine_learning/20180622/machine_learning_functions.R')
lapply(pkgs, require, character.only = T)
# registerDoParallel(cores = 20)

#---------------------------------#
#----------preprocessing----------#
#---------------------------------#
MG <- readRDS(paste0(path, 'MG_regressed_out.rds'))

mydata <- t(MG@scale.data)

# extract var genes
gene_filter_setting <- rowSums(MG@raw.data > 0) >= 100
MG_filter_genes <- row.names(MG@raw.data[gene_filter_setting,])
mydata <- mydata[,MG_filter_genes]


metadata <- data.frame(MG@ident, cell='MG')
colnames(metadata)[1]='cell_type'
metadata$cell_type <- as.character(metadata$cell_type)

metadata <- metadata[metadata$cell_type == 'Resting MG' | metadata$cell_type == "DAM",]
mydata <- mydata[rownames(mydata) %in% rownames(metadata),]
# ensure that the sample order are the same 
metadata <- metadata[rownames(mydata),]
mydata <- as.matrix(mydata)

#------------------------------------------------------#
#-----split data into training and test sets-----------#
#------------------------------------------------------#
set.seed(2018)
nall <- nrow(mydata) #total number of rows in data
ntrain <- floor(0.7 * nall) # number of rows for train,70%
ntest <- floor(0.3* nall) # number of rows for test, 30%
index <- seq(1:nall)
trainIndex <- sample(index, ntrain) #train data set
testIndex <- index[-trainIndex]
train <- mydata[trainIndex,]
test <- mydata[testIndex,]
# calculate the pre-process parameters from the dataset
preprocessParams <- preProcess(train, method=c("center", "scale"))
trainStd <- predict(preprocessParams, train)
testStd <- predict(preprocessParams, test)

mdlX <- trainStd
mdlY <- as.factor(metadata[trainIndex,][,1])
newX <- testStd
newY <- as.factor(metadata[testIndex,][,1])
learnCV <- data.frame(y = mdlY, mdlX)
testCV <-  data.frame(y = newY, newX)

# imbalance dataset to balanced dataset
set.seed(2018)
learnCV <- SMOTE(y ~., learnCV, perc.over=200, perc.under=150)



CV_result <- foreach (seed = c(2018:2022),.combine = rbind) %do% {
	# function for lasso and elastic net logistic regression based on stablity selection
	model <- enet_CV(data=learnCV, n_folds=10, model_name='glmnet.lasso', cores=10, seed)
	model
}
CV_result_lasso <- CV_result %>% group_by(cutoff) %>% summarise(Average=mean(accuracy), StDev=sd(accuracy))

# 10 folds cross validation repeated 5 times
CV_result <- foreach (seed = c(2018:2022),.combine = rbind) %do% {
	model <- enet_CV(data=learnCV, n_folds=10, model_name='glmnet.enet_0.2', cores=10, seed)
	model
}
CV_result_enet_0.2 <- CV_result %>% group_by(cutoff) %>% summarise(Average=mean(accuracy), StDev=sd(accuracy))


CV_result <- foreach (seed = c(2018:2022),.combine = rbind) %do% {
	model <- enet_CV(data=learnCV, n_folds=10, model_name='glmnet.enet_0.4', cores=10, seed)
	model
}
CV_result_enet_0.4 <- CV_result %>% group_by(cutoff) %>% summarise(Average=mean(accuracy), StDev=sd(accuracy))


CV_result <- foreach (seed = c(2018:2022),.combine = rbind) %do% {
	model <- enet_CV(data=learnCV, n_folds=10, model_name='glmnet.enet_0.6', cores=10, seed)
	model
}
CV_result_enet_0.6 <- CV_result %>% group_by(cutoff) %>% summarise(Average=mean(accuracy), StDev=sd(accuracy))


CV_result <- foreach (seed = c(2018:2022),.combine = rbind) %do% {
	model <- enet_CV(data=learnCV, n_folds=10, model_name='glmnet.enet_0.8', cores=10, seed)
	model
}
CV_result_enet_0.8 <- CV_result %>% group_by(cutoff) %>% summarise(Average=mean(accuracy), StDev=sd(accuracy))

write.csv(CV_result_lasso, 'CV_result_lasso.csv', row.names=F)

for (i in c(0.2,0.4,0.6,0.8)) {
	result_name <- paste0('CV_result_enet_', i)
	result <- get(result_name)
	write.csv(result, paste0(result_name,'.csv'), row.names=F)
}

# selected model: elastic net, alpha=0.2; cutoff 0f stability selection: 0.9
# set.seed(2018)
# options(expressions = 5e5)
# stab <- stabs::stabsel(x = learnCV[-1], y = learnCV[,1], fitfun = 'glmnet.lasso', args.fitfun = list(family = "binomial"), cutoff = 0.6, PFER=1, papply = parallel::mclapply, mc.cores = 20, sampling.type = "SS")
# stab$selected[2]

# svm_test <- function(trainX, trainY, testX, testY, features) {
# 	set.seed(2018)
# 	if (class(features) == 'integer') {
# 		svm <- train(trainX[features], trainY, method = "svmLinear",trControl = trainControl(method='none'))
# 		predictY <- predict(svm, testX[features], type = "raw")
# 		result <- confusionMatrix(predictY, testY)
# 		return(result)
# 	} 
# 	else {
# 		svm <- train(trainX, trainY, method = "svmLinear",trControl = trainControl(method='none'))
# 		predictY <- predict(svm, testX, type = "raw")
# 		result <- confusionMatrix(predictY, testY)
# 		return(result)
# 	}
# }

# svm_rfe_opt_variables <- read.csv('/share/home/duanshumin/wilsonmak/analysis/AD_microglia/GSE103334/machine_learning/20180621/SVM-RFE/optVariables.csv', stringsAsFactors = F)
# svm_rfe_opt_variables <- svm_rfe_opt_variables[2][1:36,]
# svm_rfe_opt_variables <- which(colnames(learnCV[-1]) %in% svm_rfe_opt_variables)

# wilcoxon_opt_variables <- read.csv('/share/home/duanshumin/wilsonmak/analysis/AD_microglia/GSE103334/machine_learning/20180622/wilcox_p_adj_26.csv', stringsAsFactors = F)$x
# wilcoxon_opt_variables <- which(colnames(learnCV[-1]) %in% wilcoxon_opt_variables)

# stability_opt_variables <- read.csv('/share/home/duanshumin/wilsonmak/analysis/AD_microglia/GSE103334/machine_learning/20180624/stability_selection/stability_selection_opt_genes.csv',stringsAsFactors=F)$x
# stability_opt_variables <- which(colnames(learnCV[-1]) %in% stability_opt_variables)

# # all features
# svm_test(learnCV[-1], learnCV[,1], testCV[-1], testCV$y,features='all')
# # svm-rfe
# svm_test(learnCV[-1], learnCV[,1], testCV[-1], testCV$y, features=svm_rfe_opt_variables)
# # wilcoxon
# svm_test(learnCV[-1], learnCV[,1], testCV[-1], testCV$y, features=wilcoxon_opt_variables)
# # stability selection
# svm_test(learnCV[-1], learnCV[,1], testCV[-1], testCV$y, features=stability_opt_variables)




# set.seed(2018)
# svm <- train(learnCV[-1][stab$selected], learnCV[,1], method = "svmLinear",trControl = trainControl(method='none'),allowParallel = TRUE)
# predictY <- predict(svm, testCV[-1][stab$selected], type = "raw")
# mean(testCV$y == predictY)
# confusionMatrix(predictY, testCV$y)