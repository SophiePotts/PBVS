# AICboost

aicboost <- function(X,
                     y,
                     family = gaussian(link = identity),
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
  
  # save family object to get link function
  fam <- family
  
  
  # specify loss functions according to the family object
  expr_loss = list(
    gaussian = expression((0.5 * (y - eta) ^ 2)),
    binomial = expression(-((y * log(exp(eta) / (1 + exp(eta)))) + ((1 - y) * log(1 - exp(eta) / (1 + exp(eta)))))),
    poisson = expression(-y * eta + exp(eta))
  )
  
  
  # weight matrices for boosting matrix
  # ((g'(eta))^2)/var(mu)
  # ((g'(eta))^2)/var(g(eta))
  # mu.eta is the first derivative of g() = mu wrt eta
  
  fW <- function(current_eta, fam) {
    
    res <- ((fam$mu.eta(current_eta)) ^ 2) / (fam$variance(fam$linkinv(current_eta)))
    return((as.numeric(res)))
  
  }
  
  
  ### log likelihood function
  
  loglik_fct <- function(fam, X = X, beta) {
    
    if (fam$family == "gaussian") {
      res <- sum(dnorm(y, mean = fam$linkinv(X %*% beta), sd = sd(y), log = T))
    }
    
    if (fam$family == "poisson") {
      res <- sum(dpois(y, lambda = fam$linkinv(X %*% beta), log = T))
    }
    
    if (fam$family == "binomial") {
      res <- sum(dbinom(y, size = 1, prob = fam$linkinv(X %*% beta), log = T))
    }
    return(res)
    
  }
  
  # get number of parameters K and number of observations n
  # and set iteration counter to 1
  K <- nrow(beta)
  n <- nrow(X)
  t_i <- 1

  
  # define hat matrices for each covariate k plus intercept
  H_j <- list()
  
  for (j in 1:(K - 1)) {
    
    H_j[[j]] <-  X[, c(1, j + 1)] %*% solve(t(X[, c(1, j + 1)]) %*% X[, c(1, j + 1)]) %*% t(X[, c(1, j + 1)])
  
  }
  
  
  # get derivative of the specified loss function
  d_eta <- deriv(expr_loss[[fam$family]], "eta")

  
  # build empty matrix for the beta vector of all t iterations
  beta_chain <- matrix(rep(NA, nrow(beta) * (t + 1)),
                       nrow = nrow(beta))
  
  # set first beta to starting values
  beta_chain[, 1] <- beta
  
  # define an intercept column of length n for evaluation of intercept-only model
  intercept_col <- matrix(rep(1, n))
  
  # set beta_0 (coefficient of intercept) in a way that mean(y) would result
  beta_chain[1, 1] <- fam$linkfun(solve(crossprod(intercept_col)) %*% t(intercept_col) %*% y)
  
  # initialize eta_chain having n rows (for each observation) and t+1 columns for each iteration
  # and set values of first iteration to X %*% beta_chain[,1] (eta of iteration 1)
  eta_chain <- matrix(NA, ncol = t + 1, nrow = n)
  eta_chain[, 1] <- X %*% beta_chain[, 1]
  
  
  if (fam$family == "gaussian" & fam$link == "log") {
    y <- log(y)
  }
  
  # initialize degrees of freedom data frame
  df_vec <- matrix(NA, ncol = t, nrow = (K - 1))
  
  # list with K-1 matrices (of dimension n*n) as list elements
  B_m <- (rep(list(matrix( NA, ncol = n, nrow = n)), K - 1))
  
  # list with t matrices (of dimension n*n) as list elements
  # saving each Boosting matrix (i.e. hat matrix) of the selected covariate
  # and initialize first Boosting matrix with intercept column
  B_m_k_star <- (rep(list(matrix(NA, ncol = n, nrow = n)), t))
  
  B_m_k_star[[1]] <-intercept_col %*% solve(crossprod(intercept_col)) %*% t(intercept_col)
  
  
  # initialize aic matrix
  aic <- matrix(NA, nrow = (K - 1), ncol = t)
  
  # k star vec to keep track which variables were already used
  k_star_vec <- c()
  

  # best iteration
  b_iter <- NA
  
  
  ###################### t loop - boosting iterations ###############
  
  
  # for each iteration
  while (t_i <= t) {
    # get current eta
    # (relevant for evaluating gradient of loss function)
    eta_chain[, t_i] <- eta <- (X %*% beta_chain[, t_i])
    
    # get current u
    u <- -attributes(eval(d_eta))$gradient
    
    # initialize vector to save beta proposals 
    beta_proposal <- matrix(NA, ncol = K - 1, nrow = 2)
    
    ################# k loop - covariate loop ##########
    
    # loop through each covariate:
    for (k in 1:(K - 1)) {
      # OLS estimate (base learner, beta_hat_j_(t)=beta_proposal)
      beta_proposal[, k] <- as.numeric(solve(t(X[, c(1, k + 1)]) %*% X[, c(1, k + 1)]) %*% t(X[, c(1, k + 1)]) %*% u)
      
      # prepare and update beta_try (using nu, (beta_tilde_j_(t))=beta_try)
      
      beta_try <- as.matrix(beta_chain[, t_i])
      
      beta_try[k + 1, 1] <- beta_try[k + 1, 1] + (nu * beta_proposal[2, k])
      
      beta_try[1, 1] <- beta_try[1, 1] + (nu * beta_proposal[1, k])
      
      
      # update degrees of freedom for k-th variable in iteration t_i
      
      df_vec[k, t_i] <- sum(diag(B_m_k_star[[t_i]])) + 
                        sum(t(nu *(fW(current_eta = eta_chain[, t_i],fam = fam) * H_j[[k]])) * (diag(n) - B_m_k_star[[t_i]]))
      
      
      # calculate AIC if k-th covariate would be used as update
      
      aic[k, t_i] <- -2 * loglik_fct(fam = fam,X = X,beta = beta_try) + 2 * df_vec[k, t_i]
      
      
    } ##### END k loop ####
    
    
    # AIC's for each possible covariate in iteration t_i
    find_min_k <- aic[, t_i]
    
    
    # save position of "best covariate",
    # i.e. the one minimizing the aic in iteration t_i
    k_star <- which.min(find_min_k)
    
    # trace k_star, the index of the "best covariate"
    k_star_vec[t_i] <- k_star
    
    # save the boosting matrix of the selected covariate in list of matrices at position of the current iteration
    B_m_k_star[[t_i + 1]] <-  B_m_k_star[[t_i]] + 
                              nu *(fW(current_eta = eta_chain[, t_i],fam = fam) * H_j[[k_star]]) %*% 
                              (diag(n) - B_m_k_star[[t_i]])
    
    # update beta chain:
    # for beta_k with nu*beta_k,
    # all other betas are held constant
    beta_chain[, t_i + 1] <- beta_chain[, t_i]
    
    beta_chain[k_star + 1, t_i + 1] <- beta_chain[k_star + 1, t_i] + (nu * beta_proposal[2, k_star])
    
    beta_chain[1, t_i + 1] <- beta_chain[1, t_i] + (nu * beta_proposal[1, k_star])
    
    # progress bar
    utils::setTxtProgressBar(progress_bar, t_i)
    
    ### check for early stopping/ best iteration and set "best performing iteration b_iter"
    if (t_i > 1) {
      
      krit1 = abs(min(aic[, t_i - 1]) - min(aic[, t_i])) / 
              abs(min(aic[, t_i]))
      
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
  # then, choose b_iter by min(aic) and replace it if it is smaller, than b_iter via krit
  if (is.na(b_iter)) {
    
    b_iter <- t_i
    
  }
  
  
  b_iter_min <- arrayInd(which.min(aic), dim(aic))[, 2]
  
  if (b_iter_min < b_iter) {
  
      b_iter <- b_iter_min
      
  }
  
  
  # calculate X*beta with last value of beta chain
  beta_end <- beta_chain[, ncol(beta_chain)]
  eta <- X %*% beta_end
  
  # apply the response function of the family object
  # to get predicted parameter of interest (e.g. lambda for poisson, pi for binomial)
  parm2pred <- fam$linkinv(eta)
  
  res <- list(beta_chain = beta_chain,
              eta = eta,
              parm2pred = parm2pred,
              aic = aic,
              df_vec = df_vec,
              B_m = B_m,
              B_m_k_star = B_m_k_star,
              eta_chain = eta_chain,
              k_star_vec = k_star_vec,
              b_iter = b_iter)
  
  return(res)
}
