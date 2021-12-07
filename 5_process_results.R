library(dplyr)

# para simu
res_para_500 <- read.csv(paste0('~/Repo/ph243/results/', "para_500_", "2021-12-02", '.csv'))
res_para_1000 <- read.csv(paste0('~/Repo/ph243/results/', "para_1000_", "2021-12-03", '.csv'))
res_sl_500 <- read.csv(paste0('~/Repo/ph243/results/', "sl_500_", "2021-12-04", '.csv'))
res_sl_1000 <- read.csv(paste0('~/Repo/ph243/results/', "sl_1000_", "2021-12-05", '.csv'))
res_hal_500 <- read.csv(paste0('~/Repo/ph243/results/', "hal_500_", "2021-12-04", '.csv'))
res_hal_1000 <- read.csv(paste0('~/Repo/ph243/results/', "hal_1000_", "2021-12-05", '.csv'))


res_para_500 <- cbind(simu_base = "para_500", res_para_500)
res_para_1000 <- cbind(simu_base = "para_1000", res_para_1000)
res_sl_500 <- cbind(simu_base = "sl_500", res_sl_500)
res_sl_1000 <-cbind(simu_base = "sl_1000", res_sl_1000)
res_hal_500 <- cbind(simu_base = "hal_500", res_hal_500)
res_hal_1000 <- cbind(simu_base = "hal_1000", res_hal_1000)



result <- dplyr::bind_rows(res_para_500 = res_para_500,
                           res_para_1000 = res_para_1000,
                           res_sl_500 = res_sl_500,
                           res_sl_1000 = res_sl_1000,
                           res_hal_500 = res_hal_500,
                           res_hal_1000 = res_hal_1000)


table_results_data <- result %>%
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
  dplyr::group_by(simu_base) %>%
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




## add hal2 and boot
library(dplyr)
res_para_boot_500 <- read.csv(paste0('~/Repo/ph243/results/', "para_boot_500_", "2021-12-06", '.csv'))

table_results_data <- res_para_boot_500 %>%
  dplyr::mutate(tmle_proportion_sl = psi_true <= tmle_upper_sl & psi_true >= tmle_lower_sl,
                cvtmle_proportion_sl = psi_true <= cvtmle_upper_sl & psi_true >= cvtmle_lower_sl,
                cvaiptw_proportion_sl = psi_true <= cvaiptw_upper_sl & psi_true >= cvaiptw_lower_sl,
                tmle_proportion_hal = psi_true <= tmle_upper_hal & psi_true >= tmle_lower_hal,
                tmle_proportion_hal2 = psi_true <= tmle_upper_hal2 & psi_true >= tmle_lower_hal2,
                tmle_proportion_boot = psi_true <= tmle_upper_boot & psi_true >= tmle_lower_boot,
                tmle_proportion_pl = psi_true <= tmle_upper_pl & psi_true >= tmle_lower_pl,
                cvtmle_proportion_hal = psi_true <= cvtmle_upper_hal & psi_true >= cvtmle_lower_hal,
                cvaiptw_proportion_hal = psi_true <= cvaiptw_upper_hal & psi_true >= cvaiptw_lower_hal,
                
                
                tmle_widthCI_sl = tmle_upper_sl-tmle_lower_sl,
                cvtmle_widthCI_sl = cvtmle_upper_sl-cvtmle_lower_sl,
                cvaiptw_widthCI_sl = cvaiptw_upper_sl-cvaiptw_lower_sl,
                tmle_widthCI_hal = tmle_upper_hal-tmle_lower_hal,
                tmle_widthCI_hal2 = tmle_upper_hal-tmle_lower_hal2,
                tmle_widthCI_boot = tmle_upper_boot-tmle_lower_boot,
                tmle_widthCI_pl = tmle_upper_boot-tmle_lower_pl,
                cvtmle_widthCI_hal = cvtmle_upper_hal-cvtmle_lower_hal,
                cvaiptw_widthCI_hal = cvaiptw_upper_hal-cvaiptw_lower_hal
  ) %>%
  dplyr::group_by(psi_true) %>%
  summarize(tmle_coverage_sl = mean(tmle_proportion_sl),
            cvtmle_coverage_sl = mean(cvtmle_proportion_sl),
            cvaiptw_coverage_sl = mean(cvaiptw_proportion_sl),
            tmle_coverage_hal = mean(tmle_proportion_hal),
            tmle_coverage_hal2 = mean(tmle_proportion_hal2),
            tmle_coverage_boot = mean(tmle_proportion_boot),
            tmle_coverage_pl = mean(tmle_proportion_pl),
            cvtmle_coverage_hal = mean(cvtmle_proportion_hal),
            cvaiptw_coverage_hal = mean(cvaiptw_proportion_hal),
            
            tmle_bias_sl = mean(tmle_sl) - mean(psi_true),
            cvtmle_bias_sl = mean(cvtmle_sl) - mean(psi_true),
            cvaiptw_bias_sl = mean(cvaiptw_sl) - mean(psi_true),
            tmle_bias_hal = mean(tmle_hal) - mean(psi_true),
            tmle_bias_hal2 = mean(tmle_hal2) - mean(psi_true),
            cvtmle_bias_hal = mean(cvtmle_hal) - mean(psi_true),
            cvaiptw_bias_hal = mean(cvaiptw_hal) - mean(psi_true),
            
            
            tmle_var_sl = var(tmle_sl),
            cvtmle_var_sl = var(cvtmle_sl),
            cvaiptw_var_sl = var(cvaiptw_sl),
            tmle_var_hal = var(tmle_hal),
            tmle_var_hal2 = var(tmle_hal2),
            cvtmle_var_hal = var(cvtmle_hal),
            cvaiptw_var_hal = var(cvaiptw_hal),
            
            tmle_mse_sl = tmle_bias_sl^2 + var(tmle_sl),
            cvtmle_mse_sl = cvtmle_bias_sl^2 + var(cvtmle_sl),
            cvaiptw_mse_sl = cvaiptw_bias_sl^2 + var(cvaiptw_sl),
            tmle_mse_hal = tmle_bias_hal^2 + var(tmle_hal),
            tmle_mse_hal2 = tmle_bias_hal2^2 + var(tmle_hal2),
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
            tmle_meanwidthCI_hal2 = mean(tmle_widthCI_hal2),
            tmle_meanwidthCI_boot = mean(tmle_widthCI_boot),
            tmle_meanwidthCI_pl = mean(tmle_widthCI_pl),
            cvtmle_meanwidthCI_hal = mean(cvtmle_widthCI_hal),
            cvaiptw_meanwidthCI_hal = mean(cvaiptw_widthCI_hal)
  )


## misQ undersmooth g

library(dplyr)
res_para_misQ_500 <- read.csv(paste0('~/Repo/ph243/results/', "para_misQ_500_", "2021-12-07", '.csv'))
res_para_misQ_1000 <- read.csv(paste0('~/Repo/ph243/results/', "para_misQ_1000_", "2021-12-07", '.csv'))
res_para_under_500 <- read.csv(paste0('~/Repo/ph243/results/', "para_under_500_", "2021-12-07", '.csv'))
res_para_under_1000 <- read.csv(paste0('~/Repo/ph243/results/', "para_under_1000_", "2021-12-07", '.csv'))



res_para_misQ_500 <- cbind(simu_base = "misQ_500", res_para_misQ_500)
res_para_misQ_1000 <- cbind(simu_base = "misQ_1000", res_para_misQ_1000)

res_para_under_500 <- cbind(simu_base = "under_500", res_para_under_500)
res_para_under_1000 <- cbind(simu_base = "under_1000", res_para_under_1000)


result <- dplyr::bind_rows(res_para_misQ_500 = res_para_misQ_500,
                           res_para_misQ_1000 = res_para_misQ_1000,
                           res_para_under_500 = res_para_under_500,
                           res_para_under_1000 = res_para_under_1000)

table_results_data <- result %>%
  dplyr::mutate(tmle_proportion_sl = psi_true <= tmle_upper_sl & psi_true >= tmle_lower_sl,
                cvtmle_proportion_sl = psi_true <= cvtmle_upper_sl & psi_true >= cvtmle_lower_sl,
                cvaiptw_proportion_sl = psi_true <= cvaiptw_upper_sl & psi_true >= cvaiptw_lower_sl,
                tmle_proportion_underhal = psi_true <= tmle_upper_underhal & psi_true >= tmle_lower_underhal,
                
                
                tmle_widthCI_sl = tmle_upper_sl-tmle_lower_sl,
                cvtmle_widthCI_sl = cvtmle_upper_sl-cvtmle_lower_sl,
                cvaiptw_widthCI_sl = cvaiptw_upper_sl-cvaiptw_lower_sl,
                tmle_widthCI_underhal= tmle_upper_underhal-tmle_lower_underhal
  ) %>%
  dplyr::group_by(simu_base) %>%
  summarize(tmle_coverage_sl = mean(tmle_proportion_sl),
            cvtmle_coverage_sl = mean(cvtmle_proportion_sl),
            cvaiptw_coverage_sl = mean(cvaiptw_proportion_sl),
            tmle_coverage_underhal = mean(tmle_proportion_underhal),

            tmle_bias_sl = mean(tmle_sl) - mean(psi_true),
            cvtmle_bias_sl = mean(cvtmle_sl) - mean(psi_true),
            cvaiptw_bias_sl = mean(cvaiptw_sl) - mean(psi_true),
            tmle_bias_underhal = mean(tmle_underhal) - mean(psi_true),
            
            tmle_var_sl = var(tmle_sl),
            cvtmle_var_sl = var(cvtmle_sl),
            cvaiptw_var_sl = var(cvaiptw_sl),
            tmle_var_underhal = var(tmle_underhal),
            
            tmle_mse_sl = tmle_bias_sl^2 + var(tmle_sl),
            cvtmle_mse_sl = cvtmle_bias_sl^2 + var(cvtmle_sl),
            cvaiptw_mse_sl = cvaiptw_bias_sl^2 + var(cvaiptw_sl),
            tmle_mse_underhal = tmle_bias_underhal^2 + var(tmle_underhal),
            
            
            # coverage of oracle CI
            # tmle_oracle = mean(psi_true <= tmle + 1.96*sd(tmle) & psi_true >= tmle - 1.96*sd(tmle)),
            # cvtmle_oracle = mean(psi_true <= cvtmle + 1.96*sd(cvtmle) & psi_true >= tmle - 1.96*sd(cvtmle)),
            # cvaiptw_oracle = mean(psi_true <= cvaiptw + 1.96*sd(cvaiptw) & psi_true >= cvaiptw - 1.96*sd(cvaiptw)),
            
            tmle_meanwidthCI_sl = mean(tmle_widthCI_sl),
            cvtmle_meanwidthCI_sl = mean(cvtmle_widthCI_sl),
            cvaiptw_meanwidthCI_sl = mean(cvaiptw_widthCI_sl),
            tmle_meanwidthCI_underhal = mean(tmle_widthCI_underhal)
  )



