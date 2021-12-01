
# ----------------------------------------
# Packages
# ----------------------------------------
library(data.table)
library(dplyr)
library(readr)
library(ggplot2)
library(origami)
library(sl3)
library(tmle3)
library(knitr)
library(kableExtra)
library(hal9001)
library(glmnet)
library(here)


# ----------------------------------------
# Sl set up
# ----------------------------------------
lrn_glm <- Lrnr_glm$new()
lrn_ridge <- Lrnr_glmnet$new(alpha = 0)
lrn_elastic <- Lrnr_glmnet$new(alpha = 0)
lrn_lasso <- Lrnr_glmnet$new(alpha = 1) 
lrn_spline <- Lrnr_polspline$new()
lrn_earth <- Lrnr_earth$new()


lrn_hal_s0 <- Lrnr_hal9001$new(max_degree = 2, reduce_basis = 0.05, smoothness_orders = 0)
lrn_hal_m0 <- Lrnr_hal9001$new(max_degree = 3, reduce_basis = 0.2, smoothness_orders = 0)
lrn_hal_m1 <- Lrnr_hal9001$new(max_degree = 3, reduce_basis = 0.2, smoothness_orders = 1)

discrete_sl_metalrn <- Lrnr_cv_selector$new()


grid_params <- list(max_depth = c(2, 4, 6, 8),
                    eta = c(0.01, 0.1, 0.2, 0.3),
                    nrounds = c(50, 100, 200))
grid <- expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)
params_default <- list(nthread = getOption("sl.cores.learners", 1))
xgb_learners <- apply(grid, MARGIN = 1, function(params_tune) {
  do.call(Lrnr_xgboost$new, c(params_default, as.list(params_tune)))})

sl_stack <- make_learner(Stack, 
                         lrn_glm,
                         lrn_elastic,
                         lrn_spline,
                         xgb_learners[[1]],
                         xgb_learners[[20]],
                         xgb_learners[[30]],
                         xgb_learners[[40]])

hal_stack <- make_learner(Stack, 
                          lrn_hal_s0,
                          lrn_hal_m0,
                          lrn_hal_m1)

# ----------------------------------------
# Data cleaning
# ----------------------------------------
washb_data <- fread(
  paste0(
    "https://raw.githubusercontent.com/tlverse/tlverse-data/master/",
    "wash-benefits/washb_data.csv"
  ),
  stringsAsFactors = TRUE
)

washb_data <- washb_data %>% 
              filter(tr %in% c("Control", "Nutrition", "Nutrition + WSH")) %>% 
              mutate(tr = case_when(tr == "Control" ~ 0,
                                    tr != "Control" ~ 1))

node_list <- list(
  W = c(
    "month", "aged", "sex", "momage", "momedu",
    "momheight", "hfiacat", "Nlt18", "Ncomp", "watmin",
    "elec", "floor", "walls", "roof", "asset_wardrobe",
    "asset_table", "asset_chair", "asset_khat",
    "asset_chouki", "asset_tv", "asset_refrig",
    "asset_bike", "asset_moto", "asset_sewmach",
    "asset_mobile"
  ),
  A = "tr",
  Y = "whz"
)

processed <- process_missing(washb_data, node_list)
washb_data <- processed$data
node_list <- processed$node_list

clean_df <- function(df = washb_data, n_covar = 3, outcome = "whz", trt = "tr"){
  # specify the outcome and covariates
  covars <- colnames(df)[-which(names(df) == outcome)]
  
  # create the sl3 task
  washb_task <- make_sl3_Task(
    data = df,
    covariates = covars,
    outcome = outcome
  )
  
  medforest <- Lrnr_ranger$new(
    num.trees = 200, write.forest = FALSE,
    importance = "impurity_corrected"
  )
  
  screen_rf <- Lrnr_screener_importance$new(learner = medforest, num_screen = n_covar)
  
  # screener must already be instantiated, we did this when we created screen_rf
  keepme <- c(trt)
  screen_augment_rf <- Lrnr_screener_augment$new(
    screener = screen_rf, default_covariates = keepme
  )
  
  # which covariates are selected on the full data?
  res_screened <- screen_augment_rf$train(washb_task)
  covars_selected <- setdiff(res_screened$fit_object$screener_selected, trt)  
  
  subdf <- df %>% dplyr::select(all_of(c(covars_selected, trt, outcome)))
  subdf <- as.data.frame(sapply(subdf, as.numeric))
  names(subdf)[names(subdf) == trt] <- 'A'
  names(subdf)[names(subdf) == outcome] <- 'Y'
  
  return(subdf)
}

# ----------------------------------------
# DGD helper functions
# ----------------------------------------
# helper function for simulating data via parametric models
simu_para <- function(n = 500){
  Z <- rbeta(n, 0.85, 0.85)
  W1 <- 4*Z - 2
  W2 <- rbinom(n, 1, 0.5)
  epsilon <- rnorm(n,0,0.25)
  
  g <- expit(W1 - W1*W2)
  A <- rbinom(n, 1, g)
  Y <- expit(W1 + A*W2 - W1*W2) + epsilon
  
  simu_data <- data.frame("Y" = Y, "A" = A, "W1" = W1, "W2" = W2)
  return(simu_data)
}

# helper function for simulating data via hal/sl fit from real data
simu_nonp <- function(df, covars, g_fit, Q_fit, rv, model_type, n = 500){
  # simulate W from the empirical distribution of W in the real data
  simu_data <- dplyr::sample_n(df, size = n, replace = TRUE)
  
  # generate predictions for A using g fit on original data and binomial sampling
  df_covars <- simu_data %>% select(all_of(covars)) %>% 
                mutate_if(sapply(., is.factor), as.numeric)
  
  if (model_type == "hal"){
    a_preds <- predict(g_fit, new_data = df_covars)
  }else {
    task_a_predict <- make_sl3_Task(data = simu_data, covariates = covars,
                                    outcome = "A", outcome_type = 'binomial')
    a_preds <- g_fit$predict(task_a_predict)
  }

  # simulate A
  simu_data$A <- rbinom(length(a_preds), 1, prob = a_preds)
  
  df_covars <- simu_data %>% select(all_of(c(covars, "A"))) %>% 
                  mutate_if(sapply(., is.factor), as.numeric)
  
  if (model_type == "hal"){
    y_preds <- predict(Q_fit, new_data = df_covars)
  }else {
    task_y_predict <- make_sl3_Task(data = simu_data, covariates = c(covars, "A"),
                                    outcome = "Y", outcome_type = 'continuous')
    y_preds <- Q_fit$predict(task_y_predict)
  }

  y_preds_error <- rnorm(length(y_preds), mean = 0, sd = sqrt(rv))
  y_preds <- y_preds + y_preds_error
  
  # simulate Y
  simu_data$Y <- y_preds
  
  # drop constant cols
  for (covs in setdiff(names(simu_data), c('Y', 'A'))) {
    if (length(unique(simu_data[[covs]])) == 1) {
      simu_data = simu_data %>% select(-all_of(covs))
    }
  }
  return(simu_data)
}


# -----------------------------------------------------------------------------
# Other helper functions
# -----------------------------------------------------------------------------

logit <- function(x){
  log(x/(1-x))
}

expit <- function(x){
  exp(x)/(1+exp(x))
}

