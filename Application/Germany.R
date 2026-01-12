### Clear workspace and load required libraries
rm(list = ls())
gc()
options(scipen = 900)

pacman::p_load(tidyverse, data.table, glmnet, sperrorest, knitr, data.table, mapview, tmap, viridisLite, Matrix, mboost, gamboostLSS, MASS, spData, sf, spdep, sphet, spatialreg, sperrorest)


cfe = function(Y,X,W) {
  ### (1) Run multiple regressions and recover residuals
  m0 = lm(Y ~ ., data = X[,-1])
  m1 = lm(as.vector(W %*% Y) ~ ., data = X[,-1])
  m2 = lm(as.vector(W  %*% W %*% Y) ~ ., data = X[,-1])
  
  e0 = m0$residuals
  e1 = m1$residuals
  e2 = m2$residuals
  
  ### (2) Estimate lambda
  a = t(e1) %*% e2
  b = t(e0) %*% e2 + t(e1) %*% e1
  c = t(e0) %*% e1
  
  D = b^2 - 4 * a * c
  
  lambda = as.numeric((b - sqrt(D)) / (2 * a))
  
  ### (3) Estimate beta
  beta = coef(m0) - lambda * coef(m1)
  
  ### (4) Estimate sigma
  sigma = as.numeric(t(e0 - lambda * e1) %*% (e0 - lambda * e1) / nrow(X))
  
  # u_hat = (diag(n) - lambda * W) %*% Y - as.matrix(X) %*% beta
  # sigma = as.numeric(crossprod(u_hat) / n)
  
  result = c(lambda = lambda, beta, sigma = sigma)
  
  return(result)
  
}

sarboost = function(Y, X, W, mstop, trace, type) {
  n = length(Y)
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
  cvr = cvrisk(mod, folds = cv_spat(map = krs_spdf, weights = model.weights(mod), type = type, B = 10))
  beta = coef(mod[mstop(cvr)], off2int = TRUE)
  
  sigma = mean(mod$resid()^2)
  
  result = c(lambda = lambda, beta, sigma = sigma)
  
  return(result)
}

sarlouis = function(Y, X, W, method = c("lasso", "adlasso")) {
  n = length(Y)
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
  
  # Calculate centroid coordinates for all geometries
  coords = suppressWarnings(st_coordinates(st_centroid(krs_spdf)))
  
  # Add centroid coordinates safely
  krs_spdf = krs_spdf %>%
    mutate(x = coords[, 1],
           y = coords[, 2]) %>%
    st_set_geometry(NULL) 
  
  # Run partitioning
  resamp = partition_kmeans(
    data = krs_spdf,
    coords = c("x", "y"),
    nfold = 10,
    repetition = 1
  )
  
  # Create a named vector of fold numbers indexed by test indices
  fold_vector = unlist(lapply(seq_along(resamp[[1]]), function(i) {
    setNames(rep(i, length(resamp[[1]][[i]]$test)), resamp[[1]][[i]]$test)
  }))
  
  if(method == "lasso") {
    cvmod = cv.glmnet(as.matrix(X[,-1]), Ystar, alpha = 1, family = "gaussian", foldid = fold_vector, standardize = TRUE)
    beta = coef(cvmod, s = "lambda.min")
    beta = as.matrix(beta)
    beta = beta[beta != 0,, drop = FALSE]  
    namebeta = rownames(beta)
    beta = as.numeric(beta)
    names(beta) = namebeta
  } else if(method == "adlasso") {
    init = ncol(X) - 1

    for (j in 2:(ncol(X) - 1 + 1)) { # assuming first column is intercept
      xj = X[, j]
      init[j-1] = sum(xj * Y) / sum(xj^2)  # x_j'y / x_j'x_j
    }
    
    # init = coef(m0) - lambda * coef(m1)
    # init = init[-1]
    
    cvmod = cv.glmnet(as.matrix(X[,-1]), Ystar, family = "gaussian", standardize = TRUE, alpha = 1, foldid = fold_vector, penalty.factor = 1 / abs(init)^3)
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

cv_spat = function(map, weights, type = "k-means spatial clustering", B = 10) {
  
  # Calculate centroid coordinates for all geometries
  coords = suppressWarnings(st_coordinates(st_centroid(map)))
  
  # Add centroid coordinates safely
  map = map %>%
    mutate(x = coords[, 1],
           y = coords[, 2]) %>%
    st_set_geometry(NULL) 
  
  # Run spatial partitioning
  resamp = partition_kmeans(
    data = map,
    coords = c("x", "y"),
    nfold = B,
    repetition = 1
  )
  
  # Initialize folds matrix
  folds = matrix(weights, nrow = length(weights), ncol = B)
  
  # Apply test indices
  for (fold in seq_len(B)) {
    test = resamp[[1]][[fold]]$test
    folds[test, fold] = 0
  }
  
  attr(folds, "type") = paste(B, "-fold ", type, sep = "")
  
  return(folds)
}

### Set seed
set.seed(123456789)

### Load Data (INKAR version 2021 (inkar_2021.csv) is freely available at https://www.inkar.de/ (DL-DE BY 2.0))
inkar = fread("Application/inkar_2021.csv")
krs_spdf = st_read(dsn = "Application/vg5000_ebenen_0101", layer = "VG5000_KRS")

### Filter for Kreise and Year of interest for cross sectional data
krs = inkar %>% filter(Raumbezug == "Kreise")
krs = krs %>% filter(Zeitbezug == 2019)

# Create wide format of Kreise for variables 
krs = krs %>%
  distinct(Name, Zeitbezug, Indikator, .keep_all = TRUE) %>%  
  pivot_wider(
    id_cols = c(Name, Kennziffer),                             
    names_from = Indikator,
    values_from = Wert
  )

krs = dplyr::select(krs, where(~ !anyNA(.)))
krs = krs %>% filter(!Name == "Eisenach, Stadt")
krs = krs%>%
  mutate(Kennziffer = sprintf("%05d", Kennziffer))

### Plot folds (random)
# Generate train/test splits
fold_list = map(1:10, function(i) {
  test_idx = sample(1:nrow(krs_spdf), size = nrow(krs_spdf) / 10)
  tibble(index = 1:nrow(krs_spdf),
         fold = i,
         set = factor(ifelse(index %in% test_idx, "test", "train"),
                      levels = c("train", "test")))
})

fold_df = bind_rows(fold_list)

map_folds_random = krs_spdf %>%
  mutate(index = row_number()) %>%
  right_join(fold_df %>% mutate(fold = factor(fold)), by = "index")


fold_to_plot = 2
map_one_fold = map_folds_random %>% filter(fold == fold_to_plot)

tm_shape(map_one_fold) +
  tm_polygons(
    "set",
    fill.scale = tm_scale(values = c(train = "skyblue", test = "salmon")),
    fill.legend = tm_legend(show = FALSE)
  )



### Plot folds (spatial)
# Calculate centroid coordinates for all geometries
coords = suppressWarnings(st_coordinates(st_centroid(krs_spdf)))

# Add centroid coordinates safely
map = krs_spdf %>%
  mutate(X = coords[, 1],
         Y = coords[, 2])


# Run partitioning
resamp = partition_kmeans(
  data = map %>% 
    st_set_geometry(NULL),
  coords = c("X", "Y"),
  nfold = 5,
  repetition = 1
)


# Create a named vector of fold numbers indexed by test indices
fold_vector = unlist(lapply(seq_along(resamp[[1]]), function(i) {
  setNames(rep(i, length(resamp[[1]][[i]]$test)), resamp[[1]][[i]]$test)
}))


fold_assignment = rep(NA, nrow(map))
fold_assignment[as.integer(names(fold_vector))] = fold_vector

map$fold = factor(fold_assignment)

fold_dfs = map_dfr(seq_along(resamp[[1]]), function(i) {
  test_idx = resamp[[1]][[i]]$test
  train_idx = resamp[[1]][[i]]$train
  bind_rows(
    tibble(index = test_idx, fold = i, set = "test"),
    tibble(index = train_idx, fold = i, set = "train")
  )
})

map_index = map %>% mutate(index = row_number())
map_folds = map_index %>%
  right_join(fold_dfs, by = "index")

fold_to_plot = 4

map_one_fold = map_folds %>% 
  filter(fold.y == fold_to_plot)

tm_shape(map_one_fold) +
  tm_polygons(
    "set",
    fill.scale = tm_scale(values = c(train = "skyblue", test = "salmon")),
    fill.legend = tm_legend(show = FALSE)
  )

### Generate spatial weight matrix
knn = knearneigh(st_centroid(krs_spdf), k = 10)
nb = knn2nb(knn, row.names = krs_spdf$AGS)
#nb = poly2nb(krs_spdf, queen = TRUE)
listw = nb2listw(nb, style = "W")
W = listw2mat(listw)  

### Filter dependent and independent variables of interest
X = krs[, c("Lebenserwartung",
            "Durchschnittsalter der Bevölkerung",
            "Arbeitslosenquote",
            "Beschäftigtenquote",
            "Erwerbsquote",
            "Selbständigenquote",
            "Schuldnerquote",
            "SGB II - Quote",
            "Beschäftigte mit akademischem Berufsabschluss", 
            "Eheschließungen",
            "Ehescheidungen",
            "Ausländeranteil",
            "Medianeinkommen",
            "Haushaltseinkommen",
            "Verbraucherinsolvenzverfahren",
            "Arbeitsvolumen",
            "Bodenfläche gesamt",
            "Wohnfläche",
            "Siedlungs- und Verkehrsfläche",
            "Waldfläche",
            "Erholungsfläche",
            "Wasserfläche",
            "Mietpreise",
            "Steuereinnahmen",
            "Bruttoinlandsprodukt je Einwohner",
            "Krankenhausversorgung",
            "Ärzte",
            "Pflegebedürftige",
            "Einwohnerdichte",
            "Pkw-Dichte",
            "Pendlersaldo",
            "Straßenverkehrsunfälle",
            "Getötete im Straßenverkehr"
)]

colnames(X) = c(
  "LIFE",
  "AGE",
  "UNEMPLOYMENT",
  "EMPLOYMENT",
  "PART",
  "SELF",
  "DEBT",
  "WELFARE",
  "ACADEMICS",
  "MARRIAGES",
  "DIVORCES",
  "FOREIGN",
  "MEDINC",
  "HHINC",
  "INS",
  "LABOR",
  "LAND",
  "LIVE",
  "URBAN",
  "FOREST",
  "RECR",
  "WATER",
  "RENT",
  "TAX",
  "GDP",
  "HOSP",
  "DR",
  "CARE",
  "POP",
  "CAR",
  "COM",
  "ACC",
  "TRF"
)



for (v in c("AGE", "MARRIAGES", "DIVORCES", "FOREIGN",
            "UNEMPLOYMENT", "EMPLOYMENT",
            "PART", "SELF", "DEBT", "WELFARE",
            "ACADEMICS","MEDINC", "HHINC", "INS", "LABOR",
            "LAND", "LIVE", "URBAN", "FOREST",
            "RECR", "WATER", "RENT", "TAX", "GDP",
            "HOSP", "DR", "CARE",
            "CAR", "POP", "COM", "ACC",
            "TRF"
)) {
  #X[[v]] = asinh(X[[v]])
  X[[v]] = scale(X[[v]], center = TRUE, scale = FALSE)
  #X[[v]] = scale(X[[v]], center = TRUE, scale = TRUE)
}

### Create design matrices
Y = scale(as.matrix(X[,1]), center = FALSE, scale = FALSE)
X = data.frame("(Intercept)" = rep(1, length(Y)), X[,-1])
colnames(X)[1] = "(Intercept)"


### Estimate the models
# (1) Maximum likelihood
mod_ml = lagsarlm(Y ~ ., data = X[,-1], listw = mat2listw(W, style = "W"), Durbin = FALSE)
mle = c(coef(mod_ml), mod_ml$s2)
names(mle) = c("lambda", names(X), "sigma")

# (2) GMM 
mod_tsls = stsls(Y ~ ., data = X[,-1], listw = mat2listw(W, style = "W"), zero.policy = FALSE, legacy = TRUE, W2X = FALSE, sig2n_k = FALSE)
gmm = c(coef(mod_tsls), crossprod(mod_tsls$residuals) / nrow(X))
names(gmm) = c("lambda", names(X), "sigma")

# (3) CFE 
ccfe = cfe(Y, X ,W)

# (4) LASSO 
LASSO = sarlouis(Y, X, W, method = "lasso")


# (5) ADLASSO
ADLASSO = sarlouis(Y,X,W, method = "adlasso")

# (6) GBM
GBM = sarboost(Y,X,W, mstop = 500, trace = FALSE, type = "k-means spatial clustering")

# Extract coefficients from all models
coefs_list = list(
  QML = mle,
  TSLS = gmm,
  CFE = ccfe,
  LASSO = LASSO,
  ADLASSO = ADLASSO,
  LTB = GBM
)


# All variable names used across models
all_vars = c(names(X), "lambda", "sigma")

# Create empty table (matrix)
coef_matrix = matrix(NA, nrow = length(all_vars), ncol = length(coefs_list),
                     dimnames = list(all_vars, names(coefs_list)))

# Fill the matrix
for (model in names(coefs_list)) {
  this_coef = coefs_list[[model]]
  coef_matrix[names(this_coef), model] = round(this_coef, 4)
}

# Convert to data.frame for nicer printing
coef_df = as.data.frame(coef_matrix)
coef_df = tibble::rownames_to_column(coef_df, var = "Variable")

# Custom variable ordering: lambda, non-W, W, sigma
non_w_vars = sort(grep("^W_", coef_df$Variable, invert = TRUE, value = TRUE))
non_w_vars = setdiff(non_w_vars, c("lambda", "sigma"))


ordered_vars = c("lambda", non_w_vars, "sigma")

coef_df = coef_df %>%
  dplyr::arrange(match(Variable, ordered_vars))

# Display the table nicely
kable(coef_df, align = "lcccccc", caption = "Coefficient estimates across different estimation strategies for German life expectancy")




