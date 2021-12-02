# ------------------------------------------------------------
# TMLE (SL/HAL), CV-TMLE (SL/HAL), CV-AIPTW (SL/HAL), ?HAL-MLE (HAL)
# ------------------------------------------------------------
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
              "ic_cvaiptw" = ic_cvaiptw)
         )
}



run_simu <- function(f_simu = simu_para, 
                     psi_true = 0.1153, 
                     n_sample = 500, 
                     N_round = 3){
  
  res_simu <- foreach(i = 1:N_round, .combine = 'rbind') %dorng% {
    cat('Starting ', i, 'th job.\n', sep = '')
    simu_data <- f_simu(n = n_sample)
    nodes <- list(W = setdiff(names(simu_data), c("A", "Y")),
                  A = "A",
                  Y = "Y")
    
    res_slest <- run_est(data = simu_data, 
                         node_list = nodes, 
                         gbound = 0.025, 
                         fancy_stack = sl_stack)
    
    res_halest <- run_est(data = simu_data, 
                          node_list = nodes, 
                          gbound = 0.025, 
                          fancy_stack = sl_stack)
    
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



