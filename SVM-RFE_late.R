library(parallel)
library(doParallel)
# cluster <- makeCluster(detectCores() - 1) # convention to leave 1 core for OS
cluster <- makeCluster(50)
registerDoParallel(cluster)


# parameters for the tune function, used for fitting the svm
trControl <- trainControl(method = "none")

# parameters for the RFE function
rfeControl <- rfeControl(functions = caretFuncs, method = "repeatedcv", number= 10, repeats=5, verbose = FALSE,allowParallel = TRUE )

set.seed(2018)
rf1 <- rfe(learnCV[-1], learnCV[,1], sizes = c(1:100) , rfeControl = rfeControl, trControl = trControl, method = "svmLinear",allowParallel = TRUE)

write.csv(rf1$result, 'SVM-RFE_top_features_line_chart_data.csv', row.names=F)
write.csv(rf1$optVariables[1:44], 'SVM-RFE_optVariables.csv', row.names=F,quote=F)