### Clear workspace and load required libraries
rm(list = ls())
gc()
options(scipen = 900)

pacman::p_load(tidyverse, bench, peakRAM, knitr, microbenchmark, Matrix, glmnet, mboost, gamboostLSS, MASS, spdep, spatialreg, furrr, future, future.apply, progressr)

plan(sequential)

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

oldboost = function(Y, X, W, mstop, trace, type) {
  ### (1) Maximize ll_con
  ll_con = function(lambda, Y, X, W) {
    n = length(Y)
    
    S = (Diagonal(n) - lambda * W)
    detS = determinant(S, logarithm = TRUE)$modulus
    
    beta = solve(t(X) %*% X) %*% t(X) %*% S %*% Y
    
    sigma = 1/n * t(Y) %*% t(S) %*% (Diagonal(n) - X %*% solve(t(X) %*% X) %*% t(X)) %*% S %*% Y
    
    nll = -as.numeric(-n/2 * (log(2 * pi) + 1) - n/2 * log(sigma) + detS)
    
    return(nll)
    
  }
  
  opt = optimize(
    f = ll_con,
    interval = c(-0.99, 0.99),
    Y = Y,
    X = as.matrix(X),
    W = W
  )
  
  lambda = opt$minimum
  
  ### (3) Boosting for SAR
  Ystar = as.vector((Diagonal(n) - lambda * W) %*% Y)
  
  mod = glmboost(Ystar ~ ., data = X, family = Gaussian(),
                 control = boost_control(trace = trace, mstop = mstop, nu = 0.1))
  cvr = cvrisk(mod, folds = cv(model.weights(mod), type = type))
  beta = coef(mod[mstop(cvr)], off2int = TRUE)
  
  sigma = mean(mod$resid()^2)
  
  result = c(lambda = lambda, beta, sigma = sigma)
  
  return(result)
}

oldlouis = function(Y, X, W, method = c("lasso", "adlasso")) {
  ### (1) Maximize ll_con
  ll_con = function(lambda, Y, X, W) {
    n = length(Y)
    
    S = (Diagonal(n) - lambda * W)
    detS = determinant(S, logarithm = TRUE)$modulus
    
    beta = solve(t(X) %*% X) %*% t(X) %*% S %*% Y
    
    sigma = 1/n * t(Y) %*% t(S) %*% (Diagonal(n) - X %*% solve(t(X) %*% X) %*% t(X)) %*% S %*% Y
    
    nll = -as.numeric(-n/2 * (log(2 * pi) + 1) - n/2 * log(sigma) + detS)
    
    return(nll)
    
  }
  
  opt = optimize(
    f = ll_con,
    interval = c(-0.99, 0.99),
    Y = Y,
    X = as.matrix(X),
    W = W
  )
  
  lambda = opt$minimum
  
  ### (3) Lasso/Adlasso for SAR
  Ystar = as.vector((diag(n) - lambda * W) %*% Y)
  
  if(method == "lasso") {
    cvmod = cv.glmnet(as.matrix(X[,-1]), Ystar, alpha = 1, family = "gaussian", standardize = TRUE)
    beta = coef(cvmod, s = "lambda.min")
    beta = as.matrix(beta)
    beta = beta[beta != 0,, drop = FALSE]  
    namebeta = rownames(beta)
    beta = as.numeric(beta)
    names(beta) = namebeta
  } else if(method == "adlasso") {
    init = as.vector(solve(t(as.matrix(X)) %*% as.matrix(X)) %*% t(as.matrix(X)) %*% Ystar)[-1]
    cvmod = cv.glmnet(as.matrix(X[,-1]), Ystar, family = "gaussian", standardize = TRUE, alpha = 1, penalty.factor = 1 / abs(init))
    beta = coef(cvmod, s = "lambda.min")
    beta = as.matrix(beta)
    beta = beta[beta != 0,, drop = FALSE]
    namebeta = rownames(beta)
    beta = as.numeric(beta)
    names(beta) = namebeta
  }
  
  animal = Ystar - as.matrix(X[, names(beta)]) %*% beta
  sigma = mean((animal)^2) 
  
  result = c(lambda = lambda, beta, sigma = sigma)
  
  return(result)
}

sarboost = function(Y, X, W, mstop, trace, type) {
  ### (1) Run regressions to obtain residuals
  WY = as.vector(W %*% Y)
  WWY = as.vector(W %*% WY)
  
  m0 = lm(Y ~ ., data = X[,-1])
  m1 = lm(WY ~ ., data = X[,-1])
  m2 = lm(WWY ~ ., data = X[,-1])
  
  e0 = m0$residuals
  e1 = m1$residuals
  e2 = m2$residuals
  
  ### (2) Estimate lambda
  a = t(e1) %*% e2
  b = t(e0) %*% e2 + t(e1) %*% e1
  c = t(e0) %*% e1
  D = b^2 - 4 * a * c
  
  lambda = as.numeric((b - sqrt(D)) / (2 * a))
  
  ### (3) Boosting for SAR
  Ystar = as.vector((Diagonal(n) - lambda * W) %*% Y)
  
  mod = glmboost(Ystar ~ ., data = X, family = Gaussian(),
                 control = boost_control(trace = trace, mstop = mstop, nu = 0.1))
  cvr = cvrisk(mod, folds = cv(model.weights(mod), type = type))
  beta = coef(mod[mstop(cvr)], off2int = TRUE)
  
  sigma = mean(mod$resid()^2)
  
  result = c(lambda = lambda, beta, sigma = sigma)
  
  return(result)
}

sarlouis = function(Y, X, W, method = c("lasso", "adlasso")) {
  ### (1) Run multiple regressions and recover residuals
  WY = as.vector(W %*% Y)
  WWY = as.vector(W %*% WY)
  
  m0 = lm(Y ~ ., data = X[,-1])
  m1 = lm(WY ~ ., data = X[,-1])
  m2 = lm(WWY ~ ., data = X[,-1])
  
  e0 = m0$residuals
  e1 = m1$residuals
  e2 = m2$residuals
  
  ### (2) Estimate lambda
  a = t(e1) %*% e2
  b = t(e0) %*% e2 + t(e1) %*% e1
  c = t(e0) %*% e1
  
  D = b^2 - 4 * a * c
  
  lambda = as.numeric((b - sqrt(D)) / (2 * a))
  
  ### (3) Lasso/Adlasso for SAR
  Ystar = as.vector((diag(n) - lambda * W) %*% Y)
  
  if(method == "lasso") {
    cvmod = cv.glmnet(as.matrix(X[,-1]), Ystar, alpha = 1, family = "gaussian", standardize = TRUE)
    beta = coef(cvmod, s = "lambda.min")
    beta = as.matrix(beta)
    beta = beta[beta != 0,, drop = FALSE]  
    namebeta = rownames(beta)
    beta = as.numeric(beta)
    names(beta) = namebeta
  } else if(method == "adlasso") {
    init = (coef(m0) - lambda * coef(m1))[-1]
    cvmod = cv.glmnet(as.matrix(X[,-1]), Ystar, family = "gaussian", standardize = TRUE, alpha = 1, penalty.factor = 1 / abs(init))
    beta = coef(cvmod, s = "lambda.min")
    beta = as.matrix(beta)
    beta = beta[beta != 0,, drop = FALSE]
    namebeta = rownames(beta)
    beta = as.numeric(beta)
    names(beta) = namebeta
  }
  
  animal = Ystar - as.matrix(X[, names(beta)]) %*% beta
  sigma = mean((animal)^2) 
  
  result = c(lambda = lambda, beta, sigma = sigma)
  
  return(result)
}

### Simulation Setup
lambda_t = 0.4
beta_t = c(1, 1, -0.5, -1.5, 1, rep(0,16))
names(beta_t) = c("(Intercept)", paste0("X", 1:(length(beta_t)-1)))
sigma_t = 1

p = length(beta_t) - 1
p_true = sum(beta_t[-1] != 0) 


sims = function(n, p, W, lambda_t, beta_t, sigma_t) {
  

X = matrix(runif(n * p, -2, 2),  nrow = n, ncol = p)
X = cbind(rep(1,n), X)
X = data.frame(X)
colnames(X) = names(beta_t)

eps = rnorm(n, mean = 0, sd = sigma_t)

Y = as.vector(solve(Diagonal(n) - lambda_t * W) %*% (as.matrix(X) %*% beta_t + eps))

nameVar = names(X[-1])[1:p]
trueVar = nameVar[1:4]
falseVar = nameVar[!nameVar %in% trueVar]

run = function(fun, ...) {
  env = new.env()

  mem = peakRAM::peakRAM({
    env$model = fun(...)
  })
  
  selectedVar = names(env$model)[!names(env$model) %in% c("lambda", "(Intercept)", "sigma")]
  true.positive = length(which(trueVar %in% selectedVar))
  false.positive = length(which(falseVar %in% selectedVar))
  
  TPR = true.positive / length(trueVar)
  TNR = 1 - false.positive / length(falseVar)
  
  list(
    model = env$model,  
    mem   = mem,
    TPR = TPR,
    TNR = TNR
  )
}

# Run all methods
results = list(
  QMLGBM    = run(oldboost, Y, X, W, mstop = 500, trace = FALSE, type = "kfold"),
  QMLLASSO  = run(oldlouis, Y, X, W, method = "lasso"),
  QMLADLASSO= run(oldlouis, Y, X, W, method = "adlasso"),
  GBM       = run(sarboost, Y, X, W, mstop = 500, trace = FALSE, type = "kfold"),
  LASSO     = run(sarlouis, Y, X, W, method = "lasso"),
  ADLASSO   = run(sarlouis, Y, X, W, method = "adlasso")
)

samson = do.call(rbind, lapply(names(results), function(m) {
  data.frame(
    Method = m,
    Time = results[[m]]$mem$Elapsed_Time_sec,
    RAM  = results[[m]]$mem$Peak_RAM_Used_MiB * 0.001048576,
    TPR  = results[[m]]$TPR,
    TNR =  results[[m]]$TNR
  )
}))

return(samson)

}


ripley = function(n_values, nsim) {
  set.seed(123456789)

  result = list()

  for (n in n_values) {
    cat("Running simulation study for n =", n, "\n")

    W = sparse_network(n, k = 3)

    # Set up progress bar
    pb = progressor(steps = nsim)

    results_for_n = future_lapply(1:nsim, function(i) {
      res = sims(n, p, W, lambda_t, beta_t, sigma_t)
      pb(sprintf("n=%d, replication %d", n, i))
      res
    }, future.seed = TRUE)

    result[[paste0("n=", n)]] = results_for_n
  }

  return(result)
}



n_values = c(100, 1000, 10000)
nsim = 100

results = ripley(n_values, nsim)


# Flatten and summarize
samson = map_dfr(names(results), function(n_name) {
  
  # Combine all replications for this n
  df_n = bind_rows(results[[n_name]])
  
  # Compute mean per Method
  df_summary = df_n %>%
    group_by(Method) %>%
    summarise(
      Time = median(Time),
      RAM  = median(RAM),
      TPR  = mean(TPR),
      TNR  = mean(TNR),
      .groups = "drop"
    ) %>%
    mutate(n = n_name)  # Add n as a column
  
  df_summary
})


samson <- samson %>%
  # Remove "n=" and convert to numeric
  mutate(n = as.numeric(gsub("n=", "", n))) %>%
  # Keep method order
  mutate(Method = factor(Method, levels = c("QMLGBM", "QMLLASSO", "QMLADLASSO", "GBM", "LASSO", "ADLASSO"))) %>%
  # Arrange by n and Method
  arrange(n, Method)

samson

