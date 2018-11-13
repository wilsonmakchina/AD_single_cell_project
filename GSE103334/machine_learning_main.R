path <- '/share/home/duanshumin/wilsonmak/analysis/AD_microglia/GSE103334/'

source('/share/home/duanshumin/wilsonmak/analysis/AD_microglia/machine_learning_functions.R')

pkgs <- list("dplyr","DMwR","caret", "glmnet", "doParallel", "foreach", "pROC", "stabs","parallel")
lapply(pkgs, require, character.only = T)

#load data into variable called mydata
rowData <- read.csv(paste0(path,'normalized_data.csv'),header=T, row.name=1)

mydata <- t(rowData)

metadata <- read.table(paste0(path,'metadata_remove_outlines.txt'), head=F, row.name=2, stringsAsFactors=FALSE)

# ctrl vs late AD
metadataCtrlAD <- metadata[metadata$V3 == 2 | metadata$V3 == 5,]
metadataCtrlAD$V1 = paste0('cluster_',metadataCtrlAD$V3)
# metadataCtrlAD <- metadata[metadata$V2 != '0w',]
#setting seed so we get same data split each time
mydata <- mydata[rownames(mydata) %in% rownames(metadataCtrlAD),]

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

# metadataCtrlAD[metadataCtrlAD[,4] == 3,][,4]=5
mdlX <- trainStd
mdlY <- as.factor(metadataCtrlAD[trainIndex,][,1])
newX <- testStd
newY <- as.factor(metadataCtrlAD[testIndex,][,1])
learnCV <- data.frame(y = mdlY, mdlX)
testCV <-  data.frame(y = newY, newX)

# imbalance dataset to balanced dataset
set.seed(2018)
learnCV <- SMOTE(y ~., learnCV, perc.over=150, perc.under=200)
learnCV[,1] <- as.character(learnCV[,1])


CV_result <- foreach (seed = c(2018:2022),.combine = rbind) %do% {
	# function for lasso and elastic net logistic regression based on stablity selection
	model <- enet_CV(data=learnCV, n_folds=10, model_name='glmnet.lasso', cores=36, seed)
	model
}
CV_result_lasso <- CV_result %>% group_by(cutoff) %>% summarise(Average=mean(accuracy), StDev=sd(accuracy))

# 10 folds cross validation repeated 5 times
CV_result <- foreach (seed = c(2018:2022),.combine = rbind) %do% {
	model <- enet_CV(data=learnCV, n_folds=10, model_name='glmnet.enet_0.2', cores=100, seed)
	model
}
CV_result_enet_0.2 <- CV_result %>% group_by(cutoff) %>% summarise(Average=mean(accuracy), StDev=sd(accuracy))


CV_result <- foreach (seed = c(2018:2022),.combine = rbind) %do% {
	model <- enet_CV(data=learnCV, n_folds=10, model_name='glmnet.enet_0.4', cores=36, seed)
	model
}
CV_result_enet_0.4 <- CV_result %>% group_by(cutoff) %>% summarise(Average=mean(accuracy), StDev=sd(accuracy))


CV_result <- foreach (seed = c(2018:2022),.combine = rbind) %do% {
	model <- enet_CV(data=learnCV, n_folds=10, model_name='glmnet.enet_0.6', cores=36, seed)
	model
}
CV_result_enet_0.6 <- CV_result %>% group_by(cutoff) %>% summarise(Average=mean(accuracy), StDev=sd(accuracy))


CV_result <- foreach (seed = c(2018:2022),.combine = rbind) %do% {
	model <- enet_CV(data=learnCV, n_folds=10, model_name='glmnet.enet_0.8', cores=36, seed)
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
set.seed(2018)

stab <- stabs::stabsel(x = learnCV[-1], y = learnCV[,1], fitfun = 'glmnet.lasso', args.fitfun = list(family = "binomial"), cutoff = 0.6, PFER=1, papply = parallel::mclapply, mc.cores = 35, sampling.type = "SS")


svm <- train(learnCV[-1][stab$selected], learnCV[,1], method = "svmLinear",trControl = trainControl(method='none'))
predictY <- predict(svm, testCV[-1][stab$selected], type = "raw")
mean(testCV$y == predictY)



# svm <- train(learnCV[-1], learnCV[,1], method = "svmLinear",trControl = trainControl(method='none'))
# predictY <- predict(svm, testCV[-1], type = "raw")
# mean(testCV$y == predictY)
# result <- confusionMatrix(predictY, testCV$y)
# mcc(result$table)

# svm <- train(learnCV[-1], learnCV[,1], method = "rf",trControl = trainControl(method='none'))
# predictY <- predict(svm, testCV[-1], type = "raw")
# mean(testCV$y == predictY)
# result <- confusionMatrix(predictY, testCV$y)



