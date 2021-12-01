# ------------------------------------------------------------
# TMLE (SL/HAL), CV-TMLE (SL/HAL), CV-AIPTW (HAL), HAL-MLE (HAL)
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
  psi_cvaipw <- mean(a*(y-Qbar1W_cv)/g1W_cv + Qbar1W_cv) - mean((1-a)*(y-Qbar0W_cv)/(1-g1W_cv) + Qbar0W_cv)
  ic_cvaipw <- (y - QbarAW_cv) * ((a/g1W_cv) - (1 - a)/(1 - g1W_cv)) + (Qbar1W_cv - Qbar0W_cv) - psi_cvaipw
  
  
  return(list("tmle_fit" = tmle_fit,
              "cvtmle_fit" = cvtmle_fit,
              "psi_cvaipw" = psi_cvaipw,
              "ic_cvaipw" = ic_cvaipw,
              "Qbar1W" = Qbar1W,
              "Qbar0W" = Qbar0W))
}