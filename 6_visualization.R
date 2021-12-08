library(tidyverse)
library(here)
library(ggplot2)

# para simu
res_para_500 <- read.csv(paste0(here(), "/results/para_500_", "2021-12-02", '.csv'))
res_para_1000 <- read.csv(paste0(here(), "/results/para_1000_", "2021-12-03", '.csv'))
res_hal_500 <- read.csv(paste0(here(), "/results/hal_500_", "2021-12-04", '.csv'))
res_hal_1000 <- read.csv(paste0(here(), "/results/hal_1000_", "2021-12-05", '.csv'))
res_sl_500 <- read.csv(paste0(here(), "/results/sl_500_", "2021-12-04", '.csv'))
res_sl_1000 <- read.csv(paste0(here(), "/results/sl_1000_", "2021-12-05", '.csv'))

table_results_data <- rbind(res_para_500, res_sl_500, res_hal_500,
                            res_para_1000, res_sl_1000, res_hal_1000) %>%
  dplyr::mutate(simu_set = rep(c("parametric","SL","HAL"), each = 500, 2),
                    sample_size = rep(c(500, 1000), each = 1500)) %>%
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
                cvaiptw_widthCI_hal = cvaiptw_upper_hal-cvaiptw_lower_hal,
                
                tmle_resid_sl = tmle_sl-psi_true,
                cvtmle_resid_sl = cvtmle_sl-psi_true,
                cvaiptw_resid_sl = cvaiptw_sl-psi_true,
                tmle_resid_hal = tmle_hal-psi_true,
                cvtmle_resid_hal = cvtmle_hal-psi_true,
                cvaiptw_resid_hal = cvaiptw_hal-psi_true
  )

visulize_data <- table_results_data %>% dplyr::group_by(simu_set, sample_size) %>%  
  dplyr::summarize(tmle_sl_coverage = mean(tmle_proportion_sl),
                   cvtmle_sl_coverage = mean(cvtmle_proportion_sl),
                   cvaiptw_sl_coverage = mean(cvaiptw_proportion_sl),
                   tmle_hal_coverage = mean(tmle_proportion_hal),
                   cvtmle_hal_coverage = mean(cvtmle_proportion_hal),
                   cvaiptw_hal_coverage = mean(cvaiptw_proportion_hal),
                   
                   tmle_sl_mse = mean(tmle_resid_sl^2),
                   cvtmle_sl_mse = mean(cvtmle_resid_sl^2),
                   cvaiptw_sl_mse = mean(cvaiptw_resid_sl^2),
                   tmle_hal_mse = mean(tmle_resid_hal^2),
                   cvtmle_hal_mse = mean(cvtmle_resid_hal^2),
                   cvaiptw_hal_mse = mean(cvaiptw_resid_hal^2),
                   
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
  )  %>% pivot_longer(!c(simu_set, sample_size), names_to = "name", values_to = "value")  %>%
  dplyr::mutate(method = rep(c("TMLE-SL", "CV-TMLE-SL", "CV-AIPTW-SL", "TMLE-HAL", "CV-TMLE-HAL", "CV-AIPTW-HAL"),6),
                name = rep(c("coverage","MSE", "CIwidth"), each = 6,2)) %>% 
  pivot_wider(names_from = name, values_from = value)


# coverage plot
visulize_data %>% dplyr::group_by(simu_set, sample_size) %>% 
  ggplot(aes(x = coverage, y = method)) +
  geom_point() + geom_vline(xintercept = 0.95, lty = 2) +
  geom_segment(aes(x = coverage-1.96*sqrt(coverage*(1-coverage)/500), 
                   xend = coverage+1.96*sqrt(coverage*(1-coverage)/500),
                   y = method, yend = method),
               arrow = arrow(angle=90, ends = "both", length = unit(0.05,"inches"))) +
  facet_grid(sample_size ~ simu_set, labeller = "label_both")

# CIwidth
visulize_data %>% ggplot(aes(x = CIwidth, 
                             y = method)) +
  geom_point() + facet_grid(sample_size ~ simu_set, labeller = "label_both")

# MSE
visulize_data %>% ggplot(aes(x = MSE, 
                             y = method)) +
  geom_point() + facet_grid(sample_size ~ simu_set, labeller = "label_both")








