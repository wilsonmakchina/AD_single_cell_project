enet_CV <- function(data, n_folds, model_name, cores, seed) {
    set.seed(seed)
    n_train <- nrow(data)
    folds_i <- sample(rep(1:n_folds, length.out = n_train))
    model <- get(model_name)
    registerDoParallel(cores=cores)
    # 10 folds CV
    search_all <- foreach(k = 1:n_folds, .combine = rbind) %dopar% {        
        test_i <- which(folds_i == k)
        xy <- data
        train_xy <- xy[-test_i,]
        test_xy <- xy[test_i,]
        #
        train_x <- train_xy[-1]
        train_y <- train_xy[,1]
        #
        test_x <- test_xy[-1]
        test_y <- test_xy[,1]
        # test_y <- ifelse(test_y==levels(test_y)[1],0,1)
        cutoff <- seq(0.6, 0.9, 0.05)
        #
        # stability selection
        set.seed(seed)
        stab <- stabs::stabsel(x = train_x, y = train_y, fitfun = model, args.fitfun = list(family = "binomial"), cutoff = 0.6, PFER=1, papply = parallel::mclapply, mc.cores = 35, sampling.type = "SS")
        search <- foreach(cutoff = cutoff, .combine = rbind) %dopar% {
            selected_genes <- which(stab$max > cutoff)
            # validation using linear svm
            fit <- train(x = train_x[selected_genes], y = train_y, method = "svmLinear",trControl = trainControl(method='none'))
            yhat.stabsel <- predict(fit, test_x[selected_genes], type = "raw")
            # accuracy
            accuracy = mean(test_y == yhat.stabsel)
            print(cutoff)
            data.frame(accuracy = accuracy, cutoff = cutoff, fold = k)
        }  
        print(1)
        search                           
    }
    return(search_all)
}

a_prameter <- function (a, x, y, q, type = c("conservative", "anticonservative"), ...) {
    if (!requireNamespace("glmnet", quietly = TRUE))
        stop("Package ", sQuote("glmnet"), " needed but not available")
    if (is.data.frame(x)) {
        message("Note: ", sQuote("x"), " is coerced to a model matrix without intercept")
        x <- model.matrix(~. - 1, x)
    }
    if ("lambda" %in% names(list(...)))
        stop("It is not permitted to specify the penalty parameter ",
            sQuote("lambda"), " for lasso when used with stability selection.")
    type <- match.arg(type)
    if (type == "conservative")
        fit <- suppressWarnings(glmnet::glmnet(x, y, pmax = q, alpha = a,
            ...))
    if (type == "anticonservative")
        fit <- glmnet::glmnet(x, y, dfmax = q - 1, alpha = a, ...)
    selected <- predict(fit, type = "nonzero")
    selected <- selected[[length(selected)]]
    ret <- logical(ncol(x))
    ret[selected] <- TRUE
    names(ret) <- colnames(x)
    cf <- fit$beta
    sequence <- as.matrix(cf != 0)
    return(list(selected = ret, path = sequence))
}

glmnet.enet_0.8 <- function (x, y, q, type = c("conservative", "anticonservative"), ...) {
    #
    a = 0.8
    a_prameter(a=a,x=x,y=y,q=q,type=type)    
}

glmnet.enet_0.6 <- function (x, y, q, type = c("conservative", "anticonservative"), ...) {
    #
    a = 0.6
    a_prameter(a=a,x=x,y=y,q=q,type=type)    
}

glmnet.enet_0.4 <- function (x, y, q, type = c("conservative", "anticonservative"), ...) {
    #
    a = 0.4
    a_prameter(a=a,x=x,y=y,q=q,type=type)    
}

glmnet.enet_0.2 <- function (x, y, q, type = c("conservative", "anticonservative"), ...) {
    #
    a = 0.2
    a_prameter(a=a,x=x,y=y,q=q,type=type)    

}
glmnet.enet_0 <- function (x, y, q, type = c("conservative", "anticonservative"), ...) {
    #
    a = 0
    a_prameter(a=a,x=x,y=y,q=q,type=type)    

}

mcc <- function(tab, pos = colnames(tab)[1]) {
  if(nrow(tab) != 2 | ncol(tab) != 2) stop("A 2x2 table is needed")
  neg <- colnames(tab)[colnames(tab) != pos]
  tp <- tab[pos, pos]
  tn <- tab[neg, neg]
  fp <- tab[pos,neg]
  fn <- tab[neg, pos]
  d1 <- tp + fp
  d2 <- tp + fn
  d3 <- tn + fp
  d4 <- tn + fn
  if(d1 == 0 | d2 == 0 | d3 == 0 | d4 == 0) return(0)
  ((tp * tn) - (fp * fn))/sqrt(d1*d2*d3*d4)  
}

