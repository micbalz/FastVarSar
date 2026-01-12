### Clear workspace and load required libraries
rm(list = ls())
gc()
options(scipen = 900)

pacman::p_load(tidyverse, data.table, knitr, Matrix, glmnet, mboost, MASS, spdep, spatialreg, furrr, future, future.apply, progressr)

plan(multisession, workers = parallel::detectCores() - 1)

handlers(global = TRUE)
handlers("cli")

# Generate spatial weight matrices of networks
sparse_network = function(N, k) {
  i_list = c()
  j_list = c()
  v_list = c()
  
  for (i in 1:N) {
    for (j in 1:k) {
      ahead = (i + j - 1) %% N + 1
      behind = (i - j - 1) %% N + 1
      
      i_list = c(i_list, i, i)
      j_list = c(j_list, ahead, behind)
      v_list = c(v_list, 1, 1)
    }
  }
  
  # Create sparse matrix
  W = sparseMatrix(i = i_list, j = j_list, x = v_list, dims = c(N, N))
  
  # Normalize rows to sum to 1
  row_sums = rowSums(W)
  W = Diagonal(x = 1 / row_sums) %*% W
  
  return(W)
}

llsar = function(Y, X, W, lambda, beta, sigma) {
  n = length(Y)
  
  S = diag(n) - lambda * W
  det_S = determinant(S, logarithm = TRUE)
  log_det_S = as.numeric(det_S$modulus)
  
  ll = - n / 2 * (log(2 * pi*sigma) + 1) +
    log_det_S - (1 / (2*sigma)) * t(S %*% Y - X %*% beta) %*% (S %*% Y - X %*% beta)
  
  return(-as.numeric(ll))  
}

### Simulation Setup
nsim = 100
n = 100
lambda_t = 0.2
beta_t = c(1, 4, -2.5, -3.5, 3, rep(0,196))
names(beta_t) = c("(Intercept)", paste0("X", 1:(length(beta_t)-1)))
sigma_t = 1

p = length(beta_t) - 1
p_true = sum(beta_t[-1] != 0) 


run = function(n, lambda_t, beta_t, sigma_t, k) {
  ### Generate adjacency matrices
  W = sparse_network(n, k = k)
  
  ### Generate train data
  X = matrix(runif(n * p, -2, 2),  nrow = n, ncol = p)
  X = cbind(rep(1,n), X)
  X = data.frame(X)
  colnames(X) = names(beta_t)
  
  eps = rnorm(n, mean = 0, sd = sigma_t)
  
  Y = as.vector(solve(Diagonal(n) - lambda_t * W) %*% (as.matrix(X) %*% beta_t + eps))
  
  Z = cbind("(Intercept)" = X[,1], WY = as.matrix(W %*% Y), X[,-1])
  
  ### Generate test data
  X_test = matrix(runif(n * p, -2, 2),  nrow = n, ncol = p)
  X_test = cbind(rep(1,n), X_test)
  X_test = data.frame(X_test)
  colnames(X_test) = names(beta_t)
  
  eps_test = rnorm(n, mean = 0, sd = sigma_t)
  
  Y_test = as.vector(solve(Diagonal(n) - lambda_t * W) %*% (as.matrix(X_test) %*% beta_t + eps_test))
  
  ### Estimation
  ### (1) Search for a subset
  init = numeric(p)
  
  for (j in 2:(p+2)) { # assuming first column is intercept
    xj = Z[, j]
    init[j-1] = sum(xj * Y) / sum(xj^2)  # x_j'y / x_j'x_j
  }
  
  cvmod = cv.glmnet(as.matrix(Z[,-1]), Y, family = "gaussian", standardize = TRUE, alpha = 1, penalty.factor = 1 / abs(init))
  subs = coef(cvmod, s = "lambda.min")
  subs = as.matrix(subs)
  subs = subs[subs != 0,, drop = FALSE]
  namesubs = rownames(subs)[rownames(subs) %in% colnames(X)]
  
  Xsubs = X[ ,namesubs]
  
  ### (2) Compute lambda
  # (2.1) Run multiple regressions and recover residuals
  WY = as.vector(W %*% Y)
  WWY = as.vector(W %*% WY)
  
  m0 = lm(Y ~ ., data = Xsubs[,-1])
  m1 = lm(WY ~ ., data = Xsubs[,-1])
  m2 = lm(WWY ~ ., data = Xsubs[,-1])
  
  e0 = m0$residuals
  e1 = m1$residuals
  e2 = m2$residuals
  
  # (2.2) Estimate lambda
  a = t(e1) %*% e2
  b = t(e0) %*% e2 + t(e1) %*% e1
  c = t(e0) %*% e1
  
  D = b^2 - 4 * a * c
  
  lambda = as.numeric((b - sqrt(D)) / (2 * a))
  
  ### (3) Lasso/Adlasso for SAR
  Ystar = as.vector((diag(n) - lambda * W) %*% Y)
  
  ### Variable Selection
  # (1) LASSO
  cvmod = cv.glmnet(as.matrix(X[,-1]), Ystar, alpha = 1, family = "gaussian", standardize = TRUE)
  lasso = coef(cvmod, s = "lambda.min")
  lasso = as.matrix(lasso)
  lasso = lasso[lasso != 0,, drop = FALSE]  
  namelasso = rownames(lasso)
  lasso = as.numeric(lasso)
  names(lasso) = namelasso
  
  animal = Ystar - as.matrix(X[, names(lasso)]) %*% lasso
  sigmalasso = mean((animal)^2) 
  
  lasso = c(lambda = lambda, lasso, sigma = sigmalasso)
  
  # (2) ADLASSO
  init = numeric(p)
  
  for (j in 2:(p+1)) { # assuming first column is intercept
    xj = X[, j]
    init[j-1] = sum(xj * Y) / sum(xj^2)  # x_j'y / x_j'x_j
  }
  
  
  cvmod2 = cv.glmnet(as.matrix(X[,-1]), Ystar, family = "gaussian", standardize = TRUE, alpha = 1, penalty.factor = 1 / abs(init))
  adlasso = coef(cvmod2, s = "lambda.min")
  adlasso = as.matrix(adlasso)
  adlasso = adlasso[adlasso != 0,, drop = FALSE]  
  nameadlasso = rownames(adlasso)
  adlasso = as.numeric(adlasso)
  names(adlasso) = nameadlasso
  
  innov = Ystar - as.matrix(X[, names(adlasso)]) %*% adlasso
  sigmaadlasso = mean((innov)^2) 
  
  adlasso = c(lambda = lambda, adlasso, sigma = sigmaadlasso)
  
  # (3) Boosting
  mod = glmboost(Ystar ~ ., data = X, family = Gaussian(),
                 control = boost_control(trace = FALSE, mstop = 500, nu = 0.1))
  cvr = cvrisk(mod, folds = cv(model.weights(mod), type = "kfold"))
  bst = coef(mod[mstop(cvr)], off2int = TRUE)
  
  sigmabst = mean(mod$resid()^2)
  
  bst = c(lambda = lambda, bst, sigma = sigmabst)
  
  ### Allocate models into list
  mods = list(
    LASSO = lasso,
    ADLASSO = adlasso,
    LTB = bst
  )
  
  ### Performance Criteria
  # (1) Variable Selection
  nameVar = names(X[-1])[1:p]
  trueVar = nameVar[1:p_true]
  falseVar = nameVar[!nameVar %in% trueVar]
  
  # extract selected variable names from each model
  selected = lapply(mods, \(v) setdiff(names(v), c("lambda", "(Intercept)", "sigma")))
  
  # compute metrics in one vectorized call
  metrics = t(sapply(selected, function(sel) c(
    TPR = sum(trueVar %in% sel) / length(trueVar),
    TNR = 1 - sum(falseVar %in% sel) / length(falseVar),
    FDR = sum(falseVar %in% sel) / length(sel)
  )))
  
  # (2) Estimation
  SE = sapply(mods, function(coef) {
    mts = intersect(names(coef), names(c(beta_t)))
    mts = mts[!mts %in% c("lambda", "(Intercept)", "sigma")]
    
    sum((c(beta_t)[mts] - coef[mts])^2)
  })
  
  # (3) Prediction (Trend)
  # Y_preds = Map(function(vars, coef_vec) {
  #   coef_trend = coef_vec[vars]
  #   coef_trend = coef_trend[!names(coef_trend) %in% c("lambda","(Intercept)","sigma")]
  # 
  #   X_sel = as.matrix(X_test[, names(coef_trend), drop = FALSE])
  #   as.vector(X_sel %*% coef_trend)
  # }, selected, mods)
  # 
  # RMSEP = sapply(Y_preds, function(p) sqrt(mean((Y_test - p)^2)))
  # MAEP = sapply(Y_preds, function(p) mean(abs(Y_test - p)))
  
  nll = list()
  
  for (m in names(mods)) {
    
    vec = mods[[m]]
    
    # Extract parameters
    lambda = as.numeric(vec["lambda"])
    sigma  = as.numeric(vec["sigma"])
    
    # Identify beta coefficient names (exclude lambda and sigma)
    beta_names = setdiff(names(vec), c("lambda", "sigma"))
    
    # Subset beta vector
    beta = as.numeric(vec[beta_names])
    names(beta) = beta_names
    
    # Subset X matrix to only columns used in this model
    X_sub = as.matrix(X_test[, beta_names])
    
    # Compute log-likelihood
    nll[[m]] = llsar(Y_test, X_sub, W, lambda = lambda, beta = beta, sigma = sigma)
  }
  
  nll = sapply(nll, function(x) x)
  
  
  ### Finalize the results
  results = data.frame(
    Model    = names(mods),
    TPR      = metrics[names(mods), "TPR"],
    TNR      = metrics[names(mods), "TNR"],
    MSE       = SE[names(mods)],
    NLL = nll
  )
  rownames(results) = NULL
  
  
  return(results)
}

sims = function(lambda_values, nsim) {
  set.seed(123456789)
  
  result = list()
  
  for (k in K) {
    
    cat("Running simulation study for k =", k, "\n")
    
    pb = progressor(along = 1:nsim)
    
    results_for_k = future_lapply(1:nsim, function(i) {
      
      res = run(n, lambda_t, beta_t, sigma_t, k)
      pb(sprintf("k=%.2f, replication %d", k, i))
      res
    }, future.seed = TRUE)
    
    result[[paste0("k=", k)]] = results_for_k
  }
  return(result)
}

K = c(1, 2, 3, 5, 10, 20)
results = sims(K, nsim)

agg = lapply(results, function(simlist) {
  rbindlist(simlist)                    
}) |> lapply(function(dt) {
  setDT(dt)[, lapply(.SD, mean), by = Model]   
})

agg = lapply(agg, function(dt) {
  dt[] = lapply(names(dt), function(nm) {
    col = dt[[nm]]
    if (nm == "NLL" && is.numeric(col)) {
      round(col, 0)        
    } else if (is.numeric(col)) {
      round(col, 3)       
    } else {
      col
    }
  })
  dt
})

agg
