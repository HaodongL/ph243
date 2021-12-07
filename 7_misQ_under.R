# library(here)
# source(paste0(here(), "/0_config.R"))
# source(paste0("~/Repo/ph243/1_hal_undersmooth.R"))
#


# temp 
# simu_data <- simu_para(n = 500)
# nodes <- list(W = setdiff(names(simu_data), c("A", "Y")),
#               A = "A",
#               Y = "Y")
# data = simu_data
# node_list = nodes
# gbound = 0.025
# fancy_stack = sl_stack
# 
# temp <- run_underHAL(simu_data)


# 
# temp2 <- run_est_missQ(data = simu_data, 
#                        node_list = nodes, 
#                        gbound = 0.025, 
#                        fancy_stack = sl_stack)


# use undersmoothed HAL for g model when Q model is misspecified
run_underHAL_g <- function(df){
  
  y_u <- max(df$Y)
  y_l <- min(df$Y)
  
  df = df %>% mutate(Y = (Y - y_l)/(y_u - y_l))
  df1 <- df %>%  mutate(A = 1)
  df0 <- df %>%  mutate(A = 0)
  
  #--- fit g
  outcome <- df$A
  x <- as.matrix(df %>% select(-c("A", "Y")))
  # initial hal fit
  hal_init_g <- undersmooth_init(X=x, Y=outcome, family = "binomial")
  # undersmoothing
  # length(hal_init_g$fit_init$coefs[-1]
  #        [which(hal_init_g$fit_init$coefs[-1] != 0)])
  
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
  
  g_hat <- predict(g_fit, x)
  
  # # temp 
  # task_a_train <- make_sl3_Task(data = df, covariates = setdiff(names(df), c("A","Y")),
  #                               outcome = "A", outcome_type = 'binomial')
  # 
  # lrnr_sl_a <- make_learner(Lrnr_sl,
  #                           learners = sl_stack,
  #                           outcome_type = 'binomial',
  #                           metalearner= make_learner(
  #                             Lrnr_solnp,
  #                             loss_function = loss_loglik_binomial,
  #                             learner_function = metalearner_logistic_binomial))
  # 
  # lrnr_sl_a_fit <- lrnr_sl_a$train(task_a_train)
  # g_hat <- lrnr_sl_a_fit$predict()
  
  # bound g
  g_hat[g_hat< 0.025] <- 0.025
  g_hat[g_hat> 0.975] <- 0.975
  
  #--- fit Q
  QbarAW <- rep(mean(df$Y), nrow(df))
  Qbar1W <- rep(mean(df$Y), nrow(df))
  Qbar0W <- rep(mean(df$Y), nrow(df))
  
  # emp_mean1 <- mean(df$Y[which(df$A == 1)])
  # emp_mean0 <- mean(df$Y[which(df$A == 0)])
  # 
  # QbarAW <- sapply(df$A, function(x) {
  #                         if (x == 1) {
  #                           emp_mean1
  #                         } else {
  #                             emp_mean0
  #                           }
  #                       })
  # 
  # Qbar1W <- rep(emp_mean1, nrow(df))
  # Qbar0W <- rep(emp_mean0, nrow(df))
  
  
  # update Q
  H_AW <- df$A/g_hat - (1-df$A)/(1-g_hat)
  
  suppressWarnings({
    logitUpdate = glm(df$Y ~ -1 + H_AW, 
                      offset = qlogis(QbarAW),
                      family = 'quasibinomial')
  })
  
  epsilon <- logitUpdate$coef
  
  QbarAW_star <- plogis(qlogis(QbarAW) + epsilon*H_AW)
  Qbar1W_star <- plogis(qlogis(Qbar1W) + epsilon/g_hat)
  Qbar0W_star <- plogis(qlogis(Qbar0W) - epsilon/(1-g_hat))
  
  # transform scale back
  QbarAW_star <- QbarAW_star * (y_u - y_l) + y_l
  Qbar1W_star <- Qbar1W_star * (y_u - y_l) + y_l
  Qbar0W_star <- Qbar0W_star * (y_u - y_l) + y_l
  df$Y <- df$Y * (y_u - y_l) + y_l
  
  # calc psi and ci
  psi <- mean(Qbar1W_star) - mean(Qbar0W_star)
  ic <- H_AW*(df$Y - QbarAW_star) + (Qbar1W_star - Qbar0W_star) - psi
  
  ci_up  = psi + 1.96*sqrt(var(ic)/nrow(df))
  ci_low = psi - 1.96*sqrt(var(ic)/nrow(df))
  
  return(list("tmle_under" = psi,
              "tmle_lower_under" =  ci_low,
              "tmle_upper_under" =  ci_up,
              "epsilon" = epsilon)
  )
}

# undersmooth both Q and g
run_underHAL <- function(df){
  
  y_u <- max(df$Y)
  y_l <- min(df$Y)
  
  df = df %>% mutate(Y = (Y - y_l)/(y_u - y_l))
  df1 <- df %>%  mutate(A = 1)
  df0 <- df %>%  mutate(A = 0)
  
  #--- fit g
  outcome <- df$A
  x <- as.matrix(df %>% select(-c("A", "Y")))
  # initial hal fit
  hal_init_g <- undersmooth_init(X=x, Y=outcome, family = "binomial")
  # undersmoothing
  # length(hal_init_g$fit_init$coefs[-1]
  #        [which(hal_init_g$fit_init$coefs[-1] != 0)])
  
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
  
  g_hat <- predict(g_fit, x)
  
  # bound g
  g_hat[g_hat< 0.025] <- 0.025
  g_hat[g_hat> 0.975] <- 0.975
  
  #--- fit Q
  outcome <- df$Y
  x <- as.matrix(df %>% select(-c("Y")))
  # initial hal fit  
  hal_init_Q <- undersmooth_init(X=x, Y=outcome, family = "gaussian")
  # undersmoothing
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
  
  QbarAW <- predict(Q_fit, x)
  
  x1 <- as.matrix(df1 %>% select(-c("Y")))
  x0 <- as.matrix(df0 %>% select(-c("Y")))
  Qbar1W <- predict(Q_fit, x1)
  Qbar0W <- predict(Q_fit, x0)
  
  # emp_mean1 <- mean(df$Y[which(df$A == 1)])
  # emp_mean0 <- mean(df$Y[which(df$A == 0)])
  # 
  # QbarAW <- sapply(df$A, function(x) {
  #                         if (x == 1) {
  #                           emp_mean1
  #                         } else {
  #                             emp_mean0
  #                           }
  #                       })
  # 
  # Qbar1W <- rep(emp_mean1, nrow(df))
  # Qbar0W <- rep(emp_mean0, nrow(df))
  
  
  # update Q
  H_AW <- df$A/g_hat - (1-df$A)/(1-g_hat)
  
  suppressWarnings({
    logitUpdate = glm(df$Y ~ -1 + H_AW, 
                      offset = qlogis(QbarAW),
                      family = 'quasibinomial')
  })
  
  epsilon <- logitUpdate$coef
  
  QbarAW_star <- plogis(qlogis(QbarAW) + epsilon*H_AW)
  Qbar1W_star <- plogis(qlogis(Qbar1W) + epsilon/g_hat)
  Qbar0W_star <- plogis(qlogis(Qbar0W) - epsilon/(1-g_hat))
  
  # transform scale back
  QbarAW_star <- QbarAW_star * (y_u - y_l) + y_l
  Qbar1W_star <- Qbar1W_star * (y_u - y_l) + y_l
  Qbar0W_star <- Qbar0W_star * (y_u - y_l) + y_l
  df$Y <- df$Y * (y_u - y_l) + y_l
  
  # calc psi and ci
  psi <- mean(Qbar1W_star) - mean(Qbar0W_star)
  ic <- H_AW*(df$Y - QbarAW_star) + (Qbar1W_star - Qbar0W_star) - psi
  
  ci_up  = psi + 1.96*sqrt(var(ic)/nrow(df))
  ci_low = psi - 1.96*sqrt(var(ic)/nrow(df))
  
  return(list("tmle_under" = psi,
              "tmle_lower_under" =  ci_low,
              "tmle_upper_under" =  ci_up,
              "epsilon" = epsilon)
  )
}






# run sl est with misspecified Q
run_est_missQ <- function(data, node_list, gbound = 0.025, fancy_stack = sl_stack){
  
  tmle_spec = tmle_ATE(treatment_level = 1, control_level = 0)
  tmle_task <- tmle_spec$make_tmle_task(data, node_list)
  
  lrnr_sl_a <- make_learner(Lrnr_sl,
                            learners = fancy_stack,
                            outcome_type = 'binomial',
                            metalearner= discrete_sl_metalrn)
  lrnr_sl_y <- lrn_mean
  
  learner_list <- list(Y = lrnr_sl_y, A = lrnr_sl_a)
  
  bounded_factor_list <- list(
    define_lf(LF_emp, "W"),
    define_lf(LF_fit, "A", learner = learner_list[["A"]], bound = gbound),
    define_lf(LF_fit, "Y", learner = learner_list[["Y"]], type = 'mean')
  )
  
  bounded_likelihood_def <- Likelihood$new(bounded_factor_list)
  bounded_likelihood <- bounded_likelihood_def$train(tmle_task)
  
  # 1. compute CV-TMLE
  targeted_likelihood <- Targeted_Likelihood$new(bounded_likelihood)
  tmle_params <- tmle_spec$make_params(tmle_task, targeted_likelihood)
  
  cvtmle_fit <- fit_tmle3(tmle_task, targeted_likelihood,
                          tmle_params,targeted_likelihood$updater)
  
  # 2. compute TMLE
  targeted_likelihood_no_cv <- Targeted_Likelihood$new(bounded_likelihood, 
                                                       updater = list(cvtmle = FALSE))
  tmle_params <- tmle_spec$make_params(tmle_task, targeted_likelihood_no_cv)
  
  tmle_fit <- fit_tmle3(tmle_task, targeted_likelihood_no_cv, 
                        tmle_params,targeted_likelihood_no_cv$updater)
  
  # tmle_fit$updater$epsilons
  epsilon <- unlist(tmle_fit$updater$epsilons)
  
  # extract cv Q(1,W), Q(0,W) and g(1|W) to compute cv-aiptw
  ate <- tmle_params[[1]]
  cf_task1 <- ate$cf_likelihood_treatment$cf_tasks[[1]]
  cf_task0 <- ate$cf_likelihood_control$cf_tasks[[1]]
  
  QbarAW_cv <- bounded_likelihood$get_likelihood(tmle_task, "Y", fold_number="validation")
  Qbar1W_cv <- bounded_likelihood$get_likelihood(cf_task1, "Y", fold_number="validation")
  Qbar0W_cv <- bounded_likelihood$get_likelihood(cf_task0, "Y", fold_number="validation")
  g1W_cv <- bounded_likelihood$get_likelihood(cf_task1, "A", fold_number="validation")
  
  # 3. CV-AIPTW
  a <- data$A
  y <- data$Y
  psi_cvaiptw <- mean(a*(y-Qbar1W_cv)/g1W_cv + Qbar1W_cv) - mean((1-a)*(y-Qbar0W_cv)/(1-g1W_cv) + Qbar0W_cv)
  ic_cvaiptw <- (y - QbarAW_cv) * ((a/g1W_cv) - (1 - a)/(1 - g1W_cv)) + (Qbar1W_cv - Qbar0W_cv) - psi_cvaiptw
  
  return(list("tmle_fit" = tmle_fit,
              "cvtmle_fit" = cvtmle_fit,
              "psi_cvaiptw" = psi_cvaiptw,
              "ic_cvaiptw" = ic_cvaiptw,
              "epsilon" = epsilon)
  )
}

run_est <- function(data, node_list, gbound = 0.025, fancy_stack = sl_stack){
  
  tmle_spec = tmle_ATE(treatment_level = 1, control_level = 0)
  tmle_task <- tmle_spec$make_tmle_task(data, node_list)
  
  lrnr_sl_a <- make_learner(Lrnr_sl,
                            learners = fancy_stack,
                            outcome_type = 'binomial',
                            metalearner= discrete_sl_metalrn)
  lrnr_sl_y <- make_learner(Lrnr_sl,
                            learners = fancy_stack,
                            outcome_type = 'continuous',
                            metalearner = discrete_sl_metalrn)
  
  learner_list <- list(Y = lrnr_sl_y, A = lrnr_sl_a)
  
  bounded_factor_list <- list(
    define_lf(LF_emp, "W"),
    define_lf(LF_fit, "A", learner = learner_list[["A"]], bound = gbound),
    define_lf(LF_fit, "Y", learner = learner_list[["Y"]], type = 'mean')
  )
  
  bounded_likelihood_def <- Likelihood$new(bounded_factor_list)
  bounded_likelihood <- bounded_likelihood_def$train(tmle_task)
  
  # 1. compute CV-TMLE
  targeted_likelihood <- Targeted_Likelihood$new(bounded_likelihood)
  tmle_params <- tmle_spec$make_params(tmle_task, targeted_likelihood)
  
  cvtmle_fit <- fit_tmle3(tmle_task, targeted_likelihood,
                          tmle_params,targeted_likelihood$updater)
  
  # 2. compute TMLE
  targeted_likelihood_no_cv <- Targeted_Likelihood$new(bounded_likelihood, 
                                                       updater = list(cvtmle = FALSE))
  tmle_params <- tmle_spec$make_params(tmle_task, targeted_likelihood_no_cv)
  
  tmle_fit <- fit_tmle3(tmle_task, targeted_likelihood_no_cv, 
                        tmle_params,targeted_likelihood_no_cv$updater)
  
  epsilon <- unlist(tmle_fit$updater$epsilons)
  
  # extract cv Q(1,W), Q(0,W) and g(1|W) to compute cv-aiptw
  ate <- tmle_params[[1]]
  cf_task1 <- ate$cf_likelihood_treatment$cf_tasks[[1]]
  cf_task0 <- ate$cf_likelihood_control$cf_tasks[[1]]
  
  QbarAW_cv <- bounded_likelihood$get_likelihood(tmle_task, "Y", fold_number="validation")
  Qbar1W_cv <- bounded_likelihood$get_likelihood(cf_task1, "Y", fold_number="validation")
  Qbar0W_cv <- bounded_likelihood$get_likelihood(cf_task0, "Y", fold_number="validation")
  g1W_cv <- bounded_likelihood$get_likelihood(cf_task1, "A", fold_number="validation")
  
  # 3. CV-AIPTW
  a <- data$A
  y <- data$Y
  psi_cvaiptw <- mean(a*(y-Qbar1W_cv)/g1W_cv + Qbar1W_cv) - mean((1-a)*(y-Qbar0W_cv)/(1-g1W_cv) + Qbar0W_cv)
  ic_cvaiptw <- (y - QbarAW_cv) * ((a/g1W_cv) - (1 - a)/(1 - g1W_cv)) + (Qbar1W_cv - Qbar0W_cv) - psi_cvaiptw
  
  return(list("tmle_fit" = tmle_fit,
              "cvtmle_fit" = cvtmle_fit,
              "psi_cvaiptw" = psi_cvaiptw,
              "ic_cvaiptw" = ic_cvaiptw,
              "epsilon" = epsilon)
  )
}


##

run_simu_misQ <- function(psi_true = 0.1153,
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
    
    res_slest <- run_est_missQ(data = simu_data,
                               node_list = nodes,
                               gbound = 0.025,
                               fancy_stack = sl_stack)
    
    res_halest <- run_underHAL_g(simu_data)
    
    #--- SL res
    df_res_sl <- extract_res(allres = res_slest)
    names(df_res_sl) <- paste0(names(df_res_sl), "_sl")
    
    #--- HAL res
    df_res_hal <- data.frame('tmle' = res_halest$tmle_under, 
                             'epsilon' = res_halest$epsilon,
                             'tmle_lower' = res_halest$tmle_lower_under, 
                             'tmle_upper' = res_halest$tmle_upper_under)
    
    names(df_res_hal) <- paste0(names(df_res_hal), "_underhal")
    
    df_res_all <- cbind("i" = i,
                        "psi_true" = psi_true,
                        df_res_sl,
                        df_res_hal)
    df_res_all
  }
  return(res_simu)
}

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
    
    res_slest <- run_est(data = simu_data,
                         node_list = nodes,
                         gbound = 0.025,
                         fancy_stack = sl_stack)
    
    res_halest <- run_underHAL(simu_data)
    
    #--- SL res
    df_res_sl <- extract_res(allres = res_slest)
    names(df_res_sl) <- paste0(names(df_res_sl), "_sl")
    
    #--- HAL res
    df_res_hal <- data.frame('tmle' = res_halest$tmle_under, 
                             'epsilon' = res_halest$epsilon,
                             'tmle_lower' = res_halest$tmle_lower_under, 
                             'tmle_upper' = res_halest$tmle_upper_under)
    
    names(df_res_hal) <- paste0(names(df_res_hal), "_underhal")
    
    df_res_all <- cbind("i" = i,
                        "psi_true" = psi_true,
                        df_res_sl,
                        df_res_hal)
    df_res_all
  }
  return(res_simu)
}



extract_res <- function(allres){
  # TMLE
  tmle_output <- allres$tmle_fit
  tmle <- tmle_output$summary$tmle_est
  tmle_se <- tmle_output$summary$se
  tmle_lower <- tmle_output$summary$lower
  tmle_upper <- tmle_output$summary$upper
  rm(tmle_output)
  
  # CVTMLE 
  cvtmle_output <- allres$cvtmle_fit
  cvtmle <- cvtmle_output$summary$tmle_est
  cvtmle_se <- cvtmle_output$summary$se
  cvtmle_lower <- cvtmle_output$summary$lower
  cvtmle_upper <- cvtmle_output$summary$upper
  rm(cvtmle_output)
  
  # CV-AIPTW 
  cvaiptw <- allres$psi_cvaiptw
  ic_cvaiptw <- allres$ic_cvaiptw
  
  cvaiptw_se <- sd(ic_cvaiptw)
  cvaiptw_upper <- cvaiptw + (1.96 * cvaiptw_se / sqrt(length(ic_cvaiptw)))
  cvaiptw_lower <- cvaiptw - (1.96 * cvaiptw_se / sqrt(length(ic_cvaiptw)))
  
  df <- data.frame('tmle' = tmle, 
                   'epsilon' = allres$epsilon,
                   'tmle_se' = tmle_se, 
                   'tmle_lower' = tmle_lower, 
                   'tmle_upper' = tmle_upper, 
                   'cvtmle' = cvtmle, 
                   'cvtmle_se' = cvtmle_se,
                   'cvtmle_lower' = cvtmle_lower, 
                   'cvtmle_upper' = cvtmle_upper, 
                   'cvaiptw' = cvaiptw , 
                   'cvaiptw_se' = cvaiptw_se, 
                   'cvaiptw_lower' = cvaiptw_lower, 
                   'cvaiptw_upper' = cvaiptw_upper)
  return(df)
}



