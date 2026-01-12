# Fast and Scalable Variable Selection for Spatial Autoregressive Models

This repository provides access to fast and scalable variable selection techniques for spatial autoregressive models.
It includes:  

- Helper functions  
- Monte Carlo experiments under a varying spatial autoregressive parameter and varying spatial weight matrices
- Estimates for the application settings on modeling the life expectancy in German districts.

The repository serves as a foundation for replication.  

## Technical Details  

For in-depth derivations and explanations of fast and scalable variable selection for spatial autoregressive models, refer to: 

tba.

## Example 
```
require(Matrix)
require(mboost)
require(glmnet)

set.seed(123456789)

# Simulate artificial data
n = 100
beta_t = c(1, 4, -2.5, -3.5, 3, rep(0,16))
names(beta_t) = c("(Intercept)", paste0("X", 1:(length(beta_t)-1)))
sigma_t = 1

p = length(beta_t) - 1
p_true = sum(beta_t[-1] != 0) 


# Generate spatial weight matrix
W = sparse_network(n, k = 3)

# Generate covariates and error
X = matrix(runif(n * p, -2, 2),  nrow = n, ncol = p)
X = cbind(rep(1,n), X)
X = data.frame(X)
colnames(X) = names(beta_t)
  
eps = rnorm(n, mean = 0, sd = sigma_t)
  
Y = as.vector(solve(Diagonal(n) - lambda_t * W) %*% (as.matrix(X) %*% beta_t + eps))
```

  
