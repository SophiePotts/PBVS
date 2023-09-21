# CVboost

cvboost <- function(X,
                    y,
                    E = 5, # attention: this is the number of folds (indicated by F in the article)
                    seed = 1996,
                    family = gaussian(link = "identity"),
                    beta = matrix(0, ncol = 1, nrow = ncol(X)),
                    nu = 0.1,
                    t = 20) {
  
  ###################### preparations ###############
  
  
  # showing a progressbar when using the function
  progress_bar <-
    utils::txtProgressBar(
      min = 1,
      max = t,
      initial = 1,
      style = 3)
  on.exit(close(progress_bar))
  
  
  # load dplyr
  library(dplyr)
  
  
  # create df and assign each observation to a fold by using the sample-command
  set.seed(seed)
  df <- data.frame(y = y,
                   X,
                   E = sample(
                     1:E,
                     size = nrow(X),
                     replace = T,
                     prob = rep((1 / E), E)))
  
  # split the data set in the folds
  #kfolds = list of data frames containing the observations per fold
  kfolds <- list()
  
  for (f in 1:E) {
    kfolds[[f]] <- df %>% filter(E == f) %>% select(-"E")
  }
  
  # get training and test data set for each fold
  # (i.e. bind all but fold f together to a train data set and save f-th fold in test data set)
  cv_train_list <- list()
  cv_test_list <- list()
  
  for (f in 1:E) {
    cv_train_list[[f]] <- do.call(rbind.data.frame, kfolds[-f])
    cv_test_list[[f]]  <- kfolds[[f]]
  }
  
  # separate the training and test data set into covariates and outcome (X and y), respectively
  ytest <- list()
  ytrain_f <- list()
  Xtest <- list()
  Xtrain <- list()
  
  
  for (f in 1:E) {
    ytest[[f]] <- cv_test_list[[f]][, 1]
    ytrain_f[[f]] <- cv_train_list[[f]][, 1]
    Xtest[[f]] <- cv_test_list[[f]][, -1]
    Xtrain[[f]] <- cv_train_list[[f]][, -1]
  }
  
  # save family object to get link function
  fam <- family
  
  # specify loss functions according to the family object
  expr_loss = list(
    gaussian = expression((0.5 * (ytrain - eta) ^ 2)),
    binomial = expression(-((ytrain * log(exp(eta) / (1 + exp(eta)))) + ((1 - ytrain) * log(1 - exp(eta) / (1 + exp(eta)))))),
    poisson = expression(-ytrain * (eta) + exp(eta))
  )
  
  # get derivative of the specified loss function
  d_eta <- deriv(expr_loss[[fam$family]], "eta")
  
  
  # get number of parameters K
  # and set iteration counter to 1
  K <- nrow(beta)
  t_i <- 1
  
  # build empty matrix for the beta vector of all t iterations
  beta_chain <- matrix(rep(NA, nrow(beta) * (t + 1)), nrow = nrow(beta))
  
  # set first beta to starting values
  beta_chain[, 1] <- beta
  intercept_col <- matrix(rep(1, nrow(y)))
  
  # set beta_0 (coefficient of intercept) in a way that mean(y) would result
  beta_chain[1, 1] <- fam$linkfun(solve(crossprod(intercept_col)) %*% t(intercept_col) %*% y)
  
  
  
  # create vector for beta to test on test data
  beta_test <- NULL
  
  
  # initialize vector, saving which variable is used in 
  # each iteration (k_star_vec) and the mean_loss over the folds
  k_star_vec <- c()
  mean_loss <- c()
  

  # best iteration
  b_iter <- NA
  
  
  ###################### t loop - boosting iterations ###############
  
  # for each iteration
  while (t_i <= t) {
    
    # calculate eta for each X_train fold
    eta_f <- lapply(
      X = Xtrain,
      FUN = function(Xtrain) { (as.matrix(Xtrain) %*% beta_chain[, t_i])})
    
    
    # calculate u for each training fold
    u <- list()
    for (f in 1:E) {
      eta <- eta_f[[f]]
      ytrain <- ytrain_f[[f]]
      
      u[[f]] <- -attributes(eval(d_eta))$gradient
    }
    
    
    # initialize vector to save beta proposals
    beta_proposal <- matrix(NA, ncol = K - 1, nrow = 2)
    
    # initalize vector to save sum of squared residuals
    #for each beta proposal
    find_min_k <- c()
    
    # loop through each covariate:
    # calculate beta_k via OLS on training data
    # update beta_test with this beta_k*nu
    # get SSR on test data
    
    find_min_k <- data.frame()
    beta_proposal <- list()
    
    ################# k loop - covariate loop ##########
    
    ################# f loop - folds loop ##########
    
    
    for (k in 1:(K - 1)) {
      for (f in 1:E) {
        
        # OLS estimate (base learner, beta_hat_[j,-f]_(t)=beta_proposal)
        beta_proposal[[f]] <- solve(t(as.matrix(Xtrain[[f]][, c(1, k + 1)])) %*%
                                      as.matrix(Xtrain[[f]][, c(1, k + 1)])) %*%
          t(as.matrix(Xtrain[[f]][, c(1, k + 1)])) %*% as.matrix(u[[f]])
        
        # prepare and update beta_try (using nu, (beta_tilde_[j,-f]_(t))=beta_try)
        beta_try <- as.matrix(beta_chain[, t_i])
        
        beta_try[k + 1, 1]  <- beta_try[k + 1, 1]  + nu * beta_proposal[[f]][2, 1]
        
        beta_try[1, 1]    <- beta_try[1, 1]    + nu * beta_proposal[[f]][1, 1]
        
        # use test data (ytest and Xtest) to evaluate loss function with beta_try
        # for each covariate (rows) for each fold (columns)
        find_min_k[k, f] <- sum((eval(d_eta, 
                                      envir = list(ytrain = as.matrix(ytest[[f]]), eta = as.matrix(Xtest[[f]]) %*% beta_try))))
        
      } ##### END f loop ####
    } ##### END k loop ####
    
    
    # get position of best beta_k
    # (i.e. find the minimum of the evaluated loss function over all folds (rowMean))
    k_star <- as.numeric(which.min(rowMeans(find_min_k)))
    
    # trace k_star, the index of the "best covariate"
    k_star_vec[t_i] <- (k_star)
    
    # # save mean(loss) for best covariate in iteration t_i
    mean_loss[t_i] <- as.numeric(rowMeans(find_min_k[k_star, ]))
    
    # use "best covariate" information to calculate coefficient vector by using the base learner on all data (train+test)
    u <- -attributes(eval(d_eta, envir = list(eta = X %*% beta_chain[, t_i],
                                              ytrain = y)))$gradient
    
    b <- solve(t(X[, c(1, k_star + 1)]) %*% X[, c(1, k_star + 1)]) %*% t(X[, c(1, k_star + 1)]) %*% u
    
    # then, update beta chain:
    # for beta_[k_star] with nu*beta_k,
    # all other betas are held constant
    beta_chain[, t_i + 1] <- beta_chain[, t_i]
    
    beta_chain[k_star + 1, t_i + 1] <- beta_chain[k_star + 1, t_i] + (nu * b[2, 1])
    
    beta_chain[1, t_i + 1] <- beta_chain[1, t_i] + (nu * b[1, 1])
    
    # progress bar
    utils::setTxtProgressBar(progress_bar, t_i)
    
    
    ### check for early stopping/ best iteration
    if (t_i > 1) {
      krit1 = abs(mean_loss[t_i - 1] - mean_loss[t_i]) / 
              abs(mean_loss[t_i])
      krit2 = sqrt(sum((beta_chain[, t_i] - beta_chain[, t_i + 1]) ^ 2)) / 
              sqrt(sum(beta_chain[, t_i + 1] ^ 2))
      
      # if(min(c(krit1, krit2)) < 1e-8) 
      if ((min(c(krit1, krit2)) < 1e-8) && (is.na(b_iter))) {
        
        b_iter <- t_i - 1
        
      }
    }
    
    
    # set iteration counter to +1
    t_i <- t_i + 1
    
    
  } ##### END t loop ####
  
  # if krit1 and krit2 never underwent the threshold, set best iteration to final iteration
  # then, choose b_iter by min(mean_loss) and replace it if it is smaller, than b_iter via krit
  
  if (is.na(b_iter)) {
    
    b_iter <- t_i - 1
    
  }
  
  b_iter_min <- which.min(mean_loss)
  
  if (b_iter_min < b_iter) {
    
    b_iter <- b_iter_min
    
  }
  
  # calculate X*beta with last value of beta chain
  beta_end <- beta_chain[, ncol(beta_chain)]
  eta <- X %*% beta_end
  
  # apply the response function of the family object
  # to get predicted parameter of interest
  # (e.g. lambda for poisson, pi for binomial)
  parm2pred <- fam$linkinv(eta)
  
  res <- list(
    beta_chain = beta_chain,
    eta = eta,
    parm2pred = parm2pred,
    k_star_vec = k_star_vec,
    mean_loss = mean_loss,
    b_iter = b_iter)
  
  return(res)
}
