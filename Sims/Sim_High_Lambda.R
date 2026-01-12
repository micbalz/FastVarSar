### Clear workspace and load required libraries
rm(list = ls())
gc()
options(scipen = 900)

pacman::p_load(tidyverse, knitr, Matrix, glmnet, mboost, MASS, spdep, spatialreg, furrr, future, future.apply, progressr)

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

cfe = function(Y, X, W) {
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
  
  return(lambda)
}

### Simulation Setup
nsim = 100
n = 100
beta_t = c(1, 4, -2.5, -3.5, 3, rep(0,196))
names(beta_t) = c("(Intercept)", paste0("X", 1:(length(beta_t)-1)))
sigma_t = 1

p = length(beta_t) - 1
p_true = sum(beta_t[-1] != 0) 


### Generate adjacency matrices
W = sparse_network(n, k = 3)

run = function(n, lambda_t, beta_t, sigma_t, W) {
  ### Generate covariates and error
  X = matrix(runif(n * p, -2, 2),  nrow = n, ncol = p)
  X = cbind(rep(1,n), X)
  X = data.frame(X)
  colnames(X) = names(beta_t)
  
  eps = rnorm(n, mean = 0, sd = sigma_t)
  
  Y = as.vector(solve(Diagonal(n) - lambda_t * W) %*% (as.matrix(X) %*% beta_t + eps))
  
  Z = cbind("(Intercept)" = X[,1], WY = as.matrix(W %*% Y), X[,-1])
  
  ### (0) CFE on oracle
  lam00 = cfe(Y, X[,names(beta_t)[beta_t != 0]], W)
  
  ### (2) CFE on lasso
  cvmod = cv.glmnet(as.matrix(Z[,-1]), Y, alpha = 1, family = "gaussian", standardize = TRUE)
  beta0 = coef(cvmod, s = "lambda.min")
  beta0 = as.matrix(beta0)
  beta0 = beta0[beta0 != 0,, drop = FALSE]
  namebeta = rownames(beta0)
  
  X_lasso = X[, namebeta[namebeta %in% colnames(X)]]
  
  lam1 = cfe(Y, X_lasso, W)
  
  ### (3) CFE on adlasso
  init = numeric(p)
  
  for (j in 2:(p+2)) { # assuming first column is intercept
    xj = Z[, j]
    init[j-1] = sum(xj * Y) / sum(xj^2)  # x_j'y / x_j'x_j
  }
  
  cvmod = cv.glmnet(as.matrix(Z[,-1]), Y, family = "gaussian", standardize = TRUE, alpha = 1, penalty.factor = 1 / abs(init))
  beta1 = coef(cvmod, s = "lambda.min")
  beta1 = as.matrix(beta1)
  beta1 = beta1[beta1 != 0,, drop = FALSE]
  namebeta1 = rownames(beta1)
  
  X_adlasso = X[, namebeta1[namebeta1 %in% colnames(X)]]
  
  lam2 = cfe(Y, X_adlasso, W)
  
  ### (4) CFE on boosting
  mod = glmboost(Y ~ ., data = Z[,-1], family = Gaussian(),
                 control = boost_control(trace = FALSE, mstop = 500, nu = 0.1))
  cvr = cvrisk(mod, folds = cv(model.weights(mod), type = "kfold"))
  beta2 = coef(mod[mstop(cvr)], off2int = TRUE)
  namebeta2 = names(beta2)[names(beta2) %in% colnames(X)]
  
  X_boost = X[, namebeta2]
  
  lam3 = cfe(Y, X_boost, W)
  
  
  results = data.frame(OR = lam00, CFELASSO = lam1, CFEADLASSO = lam2, LTB = lam3)
  
  
  return(results)
  
}


sims = function(lambda_values, nsim) {
  set.seed(123456789)
  
  result = list()
  
  for (lambda in lambda_values) {
    
    cat("Running simulation study for lambda =", lambda, "\n")
    
    pb = progressor(along = 1:nsim)
    
    results_for_lambda = future_lapply(1:nsim, function(i) {
      
      res = run(n, lambda, beta_t, sigma_t, W)
      pb(sprintf("lambda=%.2f, replication %d", lambda, i))
      res
    }, future.seed = TRUE)
    
    result[[paste0("lambda=", lambda)]] = results_for_lambda
  }
  return(result)
}

lambda_values = c(-0.8, -0.6, -0.4, -0.2, 0.2, 0.4, 0.6, 0.8)
results = sims(lambda_values, nsim)

df = imap_dfr(results, ~ {
  # .x is the list of results for this lambda
  # .y is the lambda name
  lam_value = as.numeric(gsub("lambda=", "", .y))
  
  # Make sure each element inside .x is a data frame
  tmp = map_dfr(.x, as.data.frame)
  
  tmp$lambda = lam_value
  tmp
}) %>%
  pivot_longer(cols = c("OR", "CFELASSO", "CFEADLASSO", "LTB"), 
               names_to = "method", values_to = "value")


summary_df = df %>%
  group_by(lambda, method) %>%
  reframe(
    Bias = mean(value, na.rm = TRUE) - unique(lambda),
    MSE  = sqrt(mean((value - unique(lambda))^2, na.rm = TRUE)),
    ESE  = sd(value, na.rm = TRUE)
  )

summary_df = summary_df %>%
  mutate(across(c(Bias, MSE, ESE), ~ round(., 3)))

eff = summary_df %>%
  group_by(lambda) %>%
  mutate(RE = MSE[method == "OR"] / MSE) %>%
  ungroup()

print(eff, n = nrow(eff))
