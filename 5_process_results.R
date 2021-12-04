library(dplyr)

# para simu
res_para_500 <- read.csv(paste0('~/Repo/ph243/results/', "para_500_", "2021-12-02", '.csv'))
res_para_1000 <- read.csv(paste0('~/Repo/ph243/results/', "para_1000_", "2021-12-03", '.csv'))

table_results_data <- res_para_1000 %>%
  dplyr::mutate(tmle_proportion_sl = psi_true <= tmle_upper_sl & psi_true >= tmle_lower_sl,
                cvtmle_proportion_sl = psi_true <= cvtmle_upper_sl & psi_true >= cvtmle_lower_sl,
                cvaiptw_proportion_sl = psi_true <= cvaiptw_upper_sl & psi_true >= cvaiptw_lower_sl,
                tmle_proportion_hal = psi_true <= tmle_upper_hal & psi_true >= tmle_lower_hal,
                cvtmle_proportion_hal = psi_true <= cvtmle_upper_hal & psi_true >= cvtmle_lower_hal,
                cvaiptw_proportion_hal = psi_true <= cvaiptw_upper_hal & psi_true >= cvaiptw_lower_hal,
                
                
                tmle_widthCI_sl = tmle_upper_sl-tmle_lower_sl,
                cvtmle_widthCI_sl = cvtmle_upper_sl-cvtmle_lower_sl,
                cvaiptw_widthCI_sl = cvaiptw_upper_sl-cvaiptw_lower_sl,
                tmle_widthCI_hal = tmle_upper_hal-tmle_lower_hal,
                cvtmle_widthCI_hal = cvtmle_upper_hal-cvtmle_lower_hal,
                cvaiptw_widthCI_hal = cvaiptw_upper_hal-cvaiptw_lower_hal
  ) %>%
  dplyr::group_by(psi_true) %>%
  summarize(tmle_coverage_sl = mean(tmle_proportion_sl),
            cvtmle_coverage_sl = mean(cvtmle_proportion_sl),
            cvaiptw_coverage_sl = mean(cvaiptw_proportion_sl),
            tmle_coverage_hal = mean(tmle_proportion_hal),
            cvtmle_coverage_hal = mean(cvtmle_proportion_hal),
            cvaiptw_coverage_hal = mean(cvaiptw_proportion_hal),
            
            tmle_bias_sl = mean(tmle_sl) - mean(psi_true),
            cvtmle_bias_sl = mean(cvtmle_sl) - mean(psi_true),
            cvaiptw_bias_sl = mean(cvaiptw_sl) - mean(psi_true),
            tmle_bias_hal = mean(tmle_hal) - mean(psi_true),
            cvtmle_bias_hal = mean(cvtmle_hal) - mean(psi_true),
            cvaiptw_bias_hal = mean(cvaiptw_hal) - mean(psi_true),
            
            
            tmle_var_sl = var(tmle_sl),
            cvtmle_var_sl = var(cvtmle_sl),
            cvaiptw_var_sl = var(cvaiptw_sl),
            tmle_var_hal = var(tmle_hal),
            cvtmle_var_hal = var(cvtmle_hal),
            cvaiptw_var_hal = var(cvaiptw_hal),
            
            tmle_mse_sl = tmle_bias_sl^2 + var(tmle_sl),
            cvtmle_mse_sl = cvtmle_bias_sl^2 + var(cvtmle_sl),
            cvaiptw_mse_sl = cvaiptw_bias_sl^2 + var(cvaiptw_sl),
            tmle_mse_hal = tmle_bias_hal^2 + var(tmle_hal),
            cvtmle_mse_hal = cvtmle_bias_hal^2 + var(cvtmle_hal),
            cvaiptw_mse_hal = cvaiptw_bias_hal^2 + var(cvaiptw_hal),
            
            
            # coverage of oracle CI
            # tmle_oracle = mean(psi_true <= tmle + 1.96*sd(tmle) & psi_true >= tmle - 1.96*sd(tmle)),
            # cvtmle_oracle = mean(psi_true <= cvtmle + 1.96*sd(cvtmle) & psi_true >= tmle - 1.96*sd(cvtmle)),
            # cvaiptw_oracle = mean(psi_true <= cvaiptw + 1.96*sd(cvaiptw) & psi_true >= cvaiptw - 1.96*sd(cvaiptw)),
            
            tmle_meanwidthCI_sl = mean(tmle_widthCI_sl),
            cvtmle_meanwidthCI_sl = mean(cvtmle_widthCI_sl),
            cvaiptw_meanwidthCI_sl = mean(cvaiptw_widthCI_sl),
            tmle_meanwidthCI_hal = mean(tmle_widthCI_hal),
            cvtmle_meanwidthCI_hal = mean(cvtmle_widthCI_hal),
            cvaiptw_meanwidthCI_hal = mean(cvaiptw_widthCI_hal)
  )










