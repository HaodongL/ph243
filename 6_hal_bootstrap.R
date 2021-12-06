# HAL-TMLE with bootstrap ci for ate

# library(here)
# source(paste0(here(), "/0_config.R"))
# library(devtools)
library(gentmle2)
# devtools::load_all("/Users/haodongli/Repo/TMLEbootstrap")


# set.seed(628957)
# library(truncnorm)
# context("ateBootstrap results should not be NA")
# simulate_data <- function(n_sim, a1, a2, b1) {
#   thresholding <- function(x, min, max) pmin(pmax(x, min), max)
# 
#   W <- truncnorm::rtruncnorm(n = n_sim, a = -10, b = 10, mean = 0, sd = 4)
#   A <- rbinom(
#     n_sim,
#     size = 1,
#     prob = thresholding(.3 + 0.1 * W * sin(a2 * W), 0.3, 0.7) +
#       rnorm(n_sim, mean = 0, sd = 0.05)
#   )
#   Y <- b1 * sin(W * a1) + A + rnorm(n_sim, 0, 1)
# 
#   X_matrix_0 <- data.frame(A, W)
#   all_df <- data.frame(Y, A, W)
# 
#   # append one value of Z
#   Z <- rep(1, n_sim)
#   X_matrix <- cbind(Z, X_matrix_0)
#   all_df <- cbind(Z, all_df)
#   return(list(Y = Y, A = A, W = data.frame(W), X_matrix = X_matrix, all_df = all_df))
# }
# ################################################################################
# # simulation
# ################################################################################
# n_sim <- 2e2
# a1 <- 1
# b1 <- 3
# a2 <- .1
# INFLATE_LAMBDA <- 1
# 
# data_sim <- simulate_data(n_sim = n_sim, a1 = a1, a2 = a2, b1 = b1)

# simu_data <-simu_para(n = 500)
# data_sim <- list(Y = simu_data$Y, 
#                  A = simu_data$A, 
#                  W = simu_data %>% select(setdiff(names(simu_data), c("A", "Y"))), 
#                  X_matrix = cbind(Z = 1, simu_data %>% select(setdiff(names(simu_data), c("Y")))), 
#                  all_df = cbind(Z = 1, simu_data))
# 
# boot_output <- ateBootstrap$new(data = data_sim, family_y = "gaussian")
# boot_output$bootstrap(1e2)
# regular_ci <- boot_output$all_boot_CI()
# 
# #
# ####
# tic()
# tune_param_fit <- ateTuneHyperparam$new(
#   data = data_sim, n_bootstrap = 1e2, family_y = "gaussian"
# )
# tune_param_fit$add_lambda(lambda_grid = 10^seq(log10(boot_output$pointTMLE$Q_fit$lambda_star), -6, length.out = 10))
# testthat::test_that("plateau tuning parameter is working", {
#   testthat::expect(!is.null(
#     tune_param_fit$select_lambda_pleateau_wald(tune_param_fit$get_lambda_df())
#   ))
# })
# toc()
# # tune_param_fit$plot_width(type_CI = "wald")
# tune_param_fit$lambdaCV <- boot_output$pointTMLE$Q_fit$lambda_star
# 
# lambda_pl <- tune_param_fit$select_lambda_pleateau_wald(df_lambda_width = tune_param_fit$df_lambda_width)
# 
# boot_output_pl <- ateBootstrap$new(data = data_sim, family_y = "gaussian", lambda1 = lambda_pl)
# boot_output_pl$bootstrap(1e2)
# 
# psi_pl <- boot_output_pl$psi_bootstrap[,1]
# 
# z_pl <- (psi_pl - mean(psi_pl)) / sd(psi_pl)
# q_025pl <- quantile(z_pl, 0.025)
# q_975pl <- quantile(z_pl, 0.975)
# 
# sigma_pl <- sqrt(n_sim)*sqrt(mean((psi_pl - boot_output$pointTMLE$Psi)^2))
# psi_lower_pl <- boot_output$pointTMLE$Psi + q_025pl*sigma_pl/sqrt(n_sim)
# psi_upper_pl<- boot_output$pointTMLE$Psi + q_975pl*sigma_pl/sqrt(n_sim)

###
# simu_data <- simu_para(n = 500)
# 
# simu_data_list <- list(Y = simu_data$Y, 
#                        A = simu_data$A, 
#                        W = simu_data %>% select(setdiff(names(simu_data), c("A", "Y"))), 
#                        X_matrix = cbind(Z = 1, simu_data %>% select(setdiff(names(simu_data), c("Y")))), 
#                        all_df = cbind(Z = 1, simu_data))
# 
# boot_output <- ateBootstrap$new(data = simu_data_list, family_y = "gaussian")
# boot_output$bootstrap(2e1)
# regular_ci <- boot_output$all_boot_CI()
# psi_hal <- boot_output$Psi
# boot_ci <- regular_ci$boot
# wald_ci <- regular_ci$wald

# rm(boot_output)


# temp <- run_boothal(simu_data)
# 
# simu_data <- simu_para(n = 100)

run_boothal <- function (simu_data = simu_data){
  
  simu_data_list <- list(Y = simu_data$Y, 
                         A = simu_data$A, 
                         W = simu_data %>% select(setdiff(names(simu_data), c("A", "Y"))), 
                         X_matrix = cbind(Z = 1, simu_data %>% select(setdiff(names(simu_data), c("Y")))), 
                         all_df = cbind(Z = 1, simu_data))
  boot_output <- ateBootstrap$new(data = simu_data_list, family_y = "gaussian")
  boot_output$bootstrap(1e2)
  regular_ci <- boot_output$all_boot_CI()
  psi_hal <- boot_output$Psi
  boot_ci <- regular_ci$boot
  wald_ci <- regular_ci$wald
  
  # plateau
  n_sim <- nrow(simu_data)
  tune_param_fit <- ateTuneHyperparam$new(
    data = simu_data_list, n_bootstrap = 1e2, family_y = "gaussian"
  )
  tune_param_fit$add_lambda(lambda_grid = 10^seq(log10(boot_output$pointTMLE$Q_fit$lambda_star), -6, length.out = 10))
  testthat::test_that("plateau tuning parameter is working", {
    testthat::expect(!is.null(
      tune_param_fit$select_lambda_pleateau_wald(tune_param_fit$get_lambda_df())
    ))
  })
  
  tune_param_fit$lambdaCV <- boot_output$pointTMLE$Q_fit$lambda_star
  
  lambda_pl <- tune_param_fit$select_lambda_pleateau_wald(df_lambda_width = tune_param_fit$df_lambda_width)
  
  boot_output_pl <- ateBootstrap$new(data = simu_data_list, family_y = "gaussian", lambda1 = lambda_pl)
  boot_output_pl$bootstrap(1e2)
  
  psi_pl <- boot_output_pl$psi_bootstrap[,1]
  
  z_pl <- (psi_pl - mean(psi_pl)) / sd(psi_pl)
  q_025pl <- quantile(z_pl, 0.025)
  q_975pl <- quantile(z_pl, 0.975)
  
  sigma_pl <- sqrt(n_sim)*sqrt(mean((psi_pl - boot_output$pointTMLE$Psi)^2))
  psi_lower_pl <- boot_output$pointTMLE$Psi + q_025pl*sigma_pl/sqrt(n_sim)
  psi_upper_pl<- boot_output$pointTMLE$Psi + q_975pl*sigma_pl/sqrt(n_sim)
  
  return(list("tmle" = psi_hal,
              "tmle_lower" =  wald_ci[1],
              "tmle_upper" =  wald_ci[2],
              "tmle_lower_b" =  boot_ci[1],
              "tmle_upper_b" =  boot_ci[2],
              "tmle_lower_pl" =  psi_lower_pl,
              "tmle_upper_pl" =  psi_upper_pl)
  )
}


# overwrite
run_simu <- function(psi_true = 0.1153, 
                     n_sample = 500, 
                     N_round = 500,
                     model_type = "para",
                     df_list = NULL){
  
  res_simu <- foreach(i = 1:N_round, .combine = 'rbind') %dorng% {
    cat('Starting ', i, 'th job.\n', sep = '')
    
    if (model_type == "para"){
      simu_data <- simu_para(n = n_sample)
    } else {
      simu_data <- df_list[[i]]
    }
    
    nodes <- list(W = setdiff(names(simu_data), c("A", "Y")),
                  A = "A",
                  Y = "Y")
    
    #--- boot HAL res
    res_boothal <- run_boothal(simu_data)
    
    res_slest <- run_est(data = simu_data, 
                         node_list = nodes, 
                         gbound = 0.025, 
                         fancy_stack = sl_stack)
    
    res_halest <- run_est(data = simu_data, 
                          node_list = nodes, 
                          gbound = 0.025, 
                          fancy_stack = hal_stack)
    
    #--- SL res
    df_res_sl <- extract_res(allres = res_slest)
    names(df_res_sl) <- paste0(names(df_res_sl), "_sl")
    
    #--- HAL res
    df_res_hal <- extract_res(allres = res_halest)
    names(df_res_hal) <- paste0(names(df_res_hal), "_hal")
    
    df_res_all <- cbind("i" = i,
                        "psi_true" = psi_true,
                        df_res_sl,
                        df_res_hal)
  
    
    df_res_all <- cbind(df_res_all, 
                        tmle_hal2 = res_boothal$tmle,
                        tmle_lower_hal2 = res_boothal$tmle_lower,
                        tmle_upper_hal2 = res_boothal$tmle_upper,
                        tmle_lower_boot = res_boothal$tmle_lower_b,
                        tmle_upper_boot = res_boothal$tmle_upper_b,
                        tmle_lower_pl = res_boothal$tmle_lower_pl,
                        tmle_upper_pl = res_boothal$tmle_upper_pl)
    df_res_all
  }
  return(res_simu)
}














