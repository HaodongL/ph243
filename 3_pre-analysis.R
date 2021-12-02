library(here)
source(paste0(here(), "/0_config.R"))
source(paste0(here(), "/1_hal_undersmooth.R"))

# -----------------------------------------
# calculate true ATE in three DGD settings 
# -----------------------------------------

#----- 1. parametric
calc_true_para <- function(n = 5*10^4){
  Z <- rbeta(n, 0.85, 0.85)
  W1 <- 4*Z - 2
  W2 <- rbinom(n, 1, 0.5)
  epsilon <- rnorm(n,0,0.25)
  
  A1 <- rep(1, n)
  Y1 <- expit(W1 + A1*W2 - W1*W2) + epsilon
  
  A0 <- rep(0, n)
  Y0 <- expit(W1 + A0*W2 - W1*W2) + epsilon
  psi_true <- mean(Y1 - Y0)
  return(psi_true)
}

set.seed(123)
true_ate_para <- calc_true_para(n = 10^5)


#----- 2. undersmoothed hal
df <- clean_df(df = washb_data, 
               n_covar = 3, 
               outcome = "whz", 
               trt = "tr")

covars <- setdiff(names(df), c("Y", "A"))


# initialize the undersmoothing procedure

# fit g (you can skip this part and directly read the saved .RDS below)
outcome <- df$A
x <- as.matrix(df %>% select(-c("A", "Y")))
# initial hal fit  
hal_init_g <- undersmooth_init(X=x, Y=outcome, family = "binomial")
# undersmoothing
length(hal_init_g$fit_init$coefs[-1]
       [which(hal_init_g$fit_init$coefs[-1] != 0)]) 

undersmooth_g <- undersmooth_hal(X = x,
                                 Y = outcome,
                                 fit_init = hal_init_g$fit_init,
                                 basis_mat = hal_init_g$basis_mat,
                                 Nlam = 20,
                                 family = "binomial")

g_fit <- fit_hal(X = x,
                 Y = outcome,
                 family = "binomial",
                 num_knots = num_knots_generator(
                    max_degree = 3,
                    smoothness_orders = 0,
                    base_num_knots_0 = max(100, ceiling(sqrt(length(outcome))))
                    ),
                 fit_control = list(
                    cv_select = FALSE,
                    n_folds = 10,
                    foldid = NULL,
                    use_min = TRUE,
                    lambda.min.ratio = 1e-4,
                    prediction_bounds = "default"
                    ),
                 lambda = undersmooth_g$lambda_under)

# fit Q
outcome <- df$Y
x <- as.matrix(df %>% select(-c("Y")))
# initial hal fit  
hal_init_Q <- undersmooth_init(X=x, Y=outcome, family = "gaussian")
# undersmoothing
all(hal_init_Q$fit_init$coefs[-1] == 0 )
undersmooth_Q <- undersmooth_hal(X = x,
                                 Y = outcome,
                                 fit_init = hal_init_Q$fit_init,
                                 basis_mat = hal_init_Q$basis_mat,
                                 Nlam = 20,
                                 family = "gaussian")

Q_fit <- fit_hal(X = x,
                 Y = outcome,
                 family = "gaussian",
                 num_knots = num_knots_generator(
                   max_degree = 3,
                   smoothness_orders = 0,
                   base_num_knots_0 = max(100, ceiling(sqrt(length(outcome))))
                 ),
                 fit_control = list(
                   cv_select = FALSE,
                   n_folds = 10,
                   foldid = NULL,
                   use_min = TRUE,
                   lambda.min.ratio = 1e-4,
                   prediction_bounds = "default"
                 ),
                 lambda = undersmooth_Q$lambda_under)

# generate a large sample with this distributions and calculate the truth
large_data <- dplyr::sample_n(df, size = 10^5, replace = TRUE)

a_preds <- predict(g_fit, new_data = large_data %>% select(-c("A","Y")))

# generate A
large_data$A <- rbinom(length(a_preds), 1, prob = a_preds)

df0 <- large_data %>% mutate(A = 0)
df1 <- large_data %>% mutate(A = 1)

y_preds <- predict(Q_fit, new_data = large_data %>% select(-"Y"))
y1_preds <- predict(Q_fit, new_data = df1 %>% select(-"Y"))
y0_preds <- predict(Q_fit, new_data = df0 %>% select(-"Y"))

# ate
true_ate_hal <- mean(y1_preds - y0_preds)
rv_hal <- sum((y_preds - large_data$Y)^2)/length(y_preds)

# save the res
res_hal <- list("g_fit_hal" = g_fit,
                "Q_fit_hal" = Q_fit,
                "true_ate_hal" = true_ate_hal,
                "rv_hal" = rv_hal)

# saveRDS(res_hal, paste0(here(), "/results/res_hal.RDS"))
# res_hal <- readRDS(paste0(here(), "/results/res_hal.RDS"))


# #----- 2.5 regular hal
# df <- clean_df(df = washb_data, 
#                n_covar = 3, 
#                outcome = "whz", 
#                trt = "tr")
# 
# covars <- setdiff(names(df), c("Y", "A"))
# 
# 
# # initialize the undersmoothing procedure
# 
# # fit g (you can skip this part and directly read the saved .RDS below)
# outcome <- df$A
# x <- as.matrix(df %>% select(-c("A", "Y")))
# # initial hal fit  
# g_fit <- fit_hal(X = x, Y = outcome, smoothness_orders = 0, family = "binomial")
# 
# # fit Q
# outcome <- df$Y
# x <- as.matrix(df %>% select(-c("Y")))
# # initial hal fit  
# Q_fit <- fit_hal(X = x, Y = outcome, smoothness_orders = 0, family = "gaussian")
# 
# # generate a large sample with this distributions and calculate the truth
# large_data <- dplyr::sample_n(df, size = 10^5, replace = TRUE)
# 
# a_preds <- predict(g_fit, new_data = large_data %>% select(-c("A","Y")))
# 
# # generate A
# large_data$A <- rbinom(length(a_preds), 1, prob = a_preds)
# 
# df0 <- large_data %>% mutate(A = 0)
# df1 <- large_data %>% mutate(A = 1)
# 
# y_preds <- predict(Q_fit, new_data = large_data %>% select(-"Y"))
# y1_preds <- predict(Q_fit, new_data = df1 %>% select(-"Y"))
# y0_preds <- predict(Q_fit, new_data = df0 %>% select(-"Y"))
# 
# # ate
# true_ate_hal <- mean(y1_preds - y0_preds)
# rv_hal <- sum((y_preds - large_data$Y)^2)/length(y_preds)
# 
# # save the res
# res_hal <- list("g_fit_hal" = g_fit,
#                 "Q_fit_hal" = Q_fit,
#                 "true_ate_hal" = true_ate_hal,
#                 "rv_hal" = rv_hal)
# 
# # saveRDS(res_hal, paste0(here(), "/results/res_hal.RDS"))
# # res_hal <- readRDS(paste0(here(), "/results/res_hal.RDS"))




#----- 3. sl
df <- clean_df(df = washb_data, 
               n_covar = 3, 
               outcome = "whz", 
               trt = "tr")

covars <- setdiff(names(df), c("Y", "A"))

# fit g
task_a_train <- make_sl3_Task(data = df, covariates = covars,
                              outcome = "A", outcome_type = 'binomial')

lrnr_sl_a <- make_learner(Lrnr_sl,
                          learners = sl_stack,
                          outcome_type = 'binomial',
                          metalearner= discrete_sl_metalrn)

lrnr_sl_a_fit <- lrnr_sl_a$train(task_a_train)

# fit Q
task_y_train <- make_sl3_Task(data = df, covariates = c(covars, "A"),
                              outcome = "Y", outcome_type = 'continuous')

lrnr_sl_y <- make_learner(Lrnr_sl,
                          learners = sl_stack,
                          outcome_type = 'continuous')

lrnr_sl_y_fit <- lrnr_sl_y$train(task_y_train)


# generate a large sample with this distributions and calculate the truth
large_data <- dplyr::sample_n(df, size = 10^5, replace = TRUE)

task_a_predict <- make_sl3_Task(data = large_data, covariates = covars,
                                outcome = "A", outcome_type = 'binomial')

a_preds <- lrnr_sl_a_fit$predict(task_a_predict)

# generate A
large_data$A <- rbinom(length(a_preds), 1, prob = a_preds)

df0 <- large_data %>% mutate(A = 0)
df1 <- large_data %>% mutate(A = 1)

task_y_predict <- make_sl3_Task(data = large_data, covariates = c(covars, "A"),
                                outcome = "Y", outcome_type = 'continuous')
y_preds <- lrnr_sl_y_fit$predict(task_y_predict)

task_y1_predict <- make_sl3_Task(data = df1, covariates = c(covars, "A"),
                                outcome = "Y", outcome_type = 'continuous')
y1_preds <- lrnr_sl_y_fit$predict(task_y1_predict)

task_y0_predict <- make_sl3_Task(data = df0, covariates = c(covars, "A"),
                                outcome = "Y", outcome_type = 'continuous')
y0_preds <- lrnr_sl_y_fit$predict(task_y0_predict)

# ate
true_ate_sl <- mean(y1_preds - y0_preds)
rv_sl <- sum((y_preds - large_data$Y)^2)/length(y_preds)

# save the res
res_sl <- list("g_fit_sl" = lrnr_sl_a_fit,
                "Q_fit_sl" = lrnr_sl_y_fit,
                "true_ate_sl" = true_ate_sl,
                "rv_sl" = rv_sl)

# saveRDS(res_sl, paste0(here(), "/results/res_sl.RDS"))
# res_sl <- readRDS(paste0(here(), "/results/res_sl.RDS"))



true_ate_sl <- calc_true_para(n = 5*10^4)


tm <- glm("Y~.", data = df)
