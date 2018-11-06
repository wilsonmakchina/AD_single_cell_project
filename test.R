path <- '/share/home/duanshumin/wilsonmak/analysis/AD_microglia/GSE98969/'

pkgs <- list("dplyr","DMwR","caret", "glmnet", "doParallel", "foreach", "pROC", "stabs","parallel","tidyr")

source("/share/home/duanshumin/wilsonmak/analysis/AD_microglia/machine_learning_functions.R")

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
metadata$cell_type <- gsub(' ','_', metadata$cell_type)
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



# CV_result <- foreach (seed = c(2018:2022),.combine = rbind) %do% {
# 	# function for lasso and elastic net logistic regression based on stablity selection
# 	model <- enet_CV(data=learnCV, n_folds=10, model_name='glmnet.lasso', cores=10, seed)
# 	model
# }
# CV_result_lasso <- CV_result %>% group_by(cutoff) %>% summarise(Average=mean(accuracy), StDev=sd(accuracy))

# 10 folds cross validation repeated 5 times
# CV_result <- foreach (seed = c(2018:2022),.combine = rbind) %do% {
# 	model <- enet_CV(data=learnCV, n_folds=10, model_name='glmnet.enet_0.2', cores=10, seed)
# 	model
# }
# CV_result_enet_0.2 <- CV_result %>% group_by(cutoff) %>% summarise(Average=mean(accuracy), StDev=sd(accuracy))


# CV_result <- foreach (seed = c(2018:2022),.combine = rbind) %do% {
# 	model <- enet_CV(data=learnCV, n_folds=10, model_name='glmnet.enet_0.4', cores=10, seed)
# 	model
# }
# CV_result_enet_0.4 <- CV_result %>% group_by(cutoff) %>% summarise(Average=mean(accuracy), StDev=sd(accuracy))


# CV_result <- foreach (seed = c(2018:2022),.combine = rbind) %do% {
# 	model <- enet_CV(data=learnCV, n_folds=10, model_name='glmnet.enet_0.6', cores=10, seed)
# 	model
# }
# CV_result_enet_0.6 <- CV_result %>% group_by(cutoff) %>% summarise(Average=mean(accuracy), StDev=sd(accuracy))


# CV_result <- foreach (seed = c(2018:2022),.combine = rbind) %do% {
# 	model <- enet_CV(data=learnCV, n_folds=10, model_name='glmnet.enet_0.8', cores=10, seed)
# 	model
# }
# CV_result_enet_0.8 <- CV_result %>% group_by(cutoff) %>% summarise(Average=mean(accuracy), StDev=sd(accuracy))

# write.csv(CV_result_lasso, 'CV_result_lasso.csv', row.names=F)

# for (i in c(0.2,0.4,0.6,0.8)) {
# 	result_name <- paste0('CV_result_enet_', i)
# 	result <- get(result_name)
# 	write.csv(result, paste0(result_name,'.csv'), row.names=F)
# }

# selected model: elastic net, alpha=0.4; cutoff 0f stability selection: 0.9
set.seed(2018)
options(expressions = 5e5)
stab <- stabs::stabsel(x = learnCV[-1], y = learnCV[,1], fitfun = 'glmnet.enet_0.4', args.fitfun = list(family = "binomial"), cutoff = 0.9, PFER=1, papply = parallel::mclapply, mc.cores = 35, sampling.type = "SS")



randomly_select_genes <- sample(colnames(learnCV), 48)
randomly_select_genes <- which(colnames(learnCV[-1]) %in% randomly_select_genes)

acc <- NULL
ROC <- NULL
selected_genes <- NULL
selected_genes[[1]] <- c(1:ncol(testCV[-1]))
selected_genes[[2]] <- stab$selected
selected_genes[[3]] <- randomly_select_genes
 
control <- trainControl(method = 'none', classProbs = TRUE, summaryFunction = twoClassSummary)

for (i in 1:3) {
	svm <- train(learnCV[-1][selected_genes[[i]]], learnCV[,1], method = "svmLinear",trControl = control, metric='ROC')
	predictY <- predict(svm, testCV[-1][selected_genes[[i]]], type = "raw")
	acc[i] <- mean(testCV$y == predictY)
	svm.probs <- predict(svm, testCV[-1][selected_genes[[i]]], type = "prob")
	ROC[[i]] <- roc(response = testCV$y, predictor = svm.probs$DAM, levels = levels(testCV$y))
}



p <- ggroc(list(all_genes=ROC[[1]], ML_selected_genes=ROC[[2]], random_selected_genes=ROC[[3]])) + 
	scale_fill_brewer(palette = "Set1")

pdf('roc_curve.pdf')
print(p)
dev.off()






# save.image('GSE98969_stabs_selected_genes.Rdata')


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





