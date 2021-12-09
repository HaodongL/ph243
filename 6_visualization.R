# data tidy function
tidy_visual <- function(res) {
  require(dplyr)
  temp <- res %>% 
    mutate(tmle_proportion_sl = psi_true <= tmle_upper_sl & psi_true >= tmle_lower_sl,
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
    ) %>% summarise(tmle_coverage_sl = mean(tmle_proportion_sl),
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
                    cvaiptw_meanwidthCI_hal = mean(cvaiptw_widthCI_hal),
                    
                    tmle_Q1widthCI_sl = quantile(tmle_widthCI_sl,0.25),
                    cvtmle_Q1widthCI_sl = quantile(cvtmle_widthCI_sl,0.25),
                    cvaiptw_Q1widthCI_sl = quantile(cvaiptw_widthCI_sl,0.25),
                    tmle_Q1widthCI_hal = quantile(tmle_widthCI_hal,0.25),
                    cvtmle_Q1widthCI_hal = quantile(cvtmle_widthCI_hal,0.25),
                    cvaiptw_Q1widthCI_hal = quantile(cvaiptw_widthCI_hal,0.25),
                    
                    tmle_Q3widthCI_sl = quantile(tmle_widthCI_sl,0.75),
                    cvtmle_Q3widthCI_sl = quantile(cvtmle_widthCI_sl,0.75),
                    cvaiptw_Q3widthCI_sl = quantile(cvaiptw_widthCI_sl,0.75),
                    tmle_Q3widthCI_hal = quantile(tmle_widthCI_hal,0.75),
                    cvtmle_Q3widthCI_hal = quantile(cvtmle_widthCI_hal,0.75),
                    cvaiptw_Q3widthCI_hal = quantile(cvaiptw_widthCI_hal,0.75),
                    
                    psi_true = mean(psi_true)
      )
  temp_p1 <- temp %>% gather("method", "coverage", tmle_coverage_sl:cvaiptw_coverage_hal) %>%
    select(coverage)
  temp_p2 <- temp %>% gather("method", "MSE", tmle_mse_sl:cvaiptw_mse_hal) %>%
    select(MSE)
  temp_p3 <- temp %>% gather("method", "CIwidth", tmle_meanwidthCI_sl:cvaiptw_meanwidthCI_hal) %>%
    select(CIwidth)
  temp_p4 <- temp %>% gather("method", "bias", tmle_bias_sl:cvaiptw_bias_hal) %>%
    select(bias)
  temp_p5 <- temp %>% gather("method", "Q1CIwidth", tmle_Q1widthCI_sl:cvaiptw_Q1widthCI_hal) %>%
    select(Q1CIwidth)
  temp_p6 <- temp %>% gather("method", "Q3CIwidth", tmle_Q3widthCI_sl:cvaiptw_Q3widthCI_hal) %>%
    select(Q3CIwidth)
  
  temp_long <- cbind(
    data.frame(psi_true = rep(temp$psi_true, 6), 
               method = c("TMLE-SL", "CV-TMLE-SL", "CV-AIPTW-SL", 
                          "TMLE-HAL", "CV-TMLE-HAL", "CV-AIPTW-HAL")),
    temp_p1, temp_p2, temp_p3, temp_p4, temp_p5, temp_p6) %>% 
    mutate(relativeMSE = MSE/MSE[1], 
           relativeCIwidth = CIwidth/CIwidth[1],
           relativeBias2 = bias^2/MSE[1])
  return(temp_long)
}

library(tidyverse)
library(here)
# simu
res_para_500 <- read.csv(paste0(here(), "/results/para_500_", "2021-12-02", '.csv'))
res_para_1000 <- read.csv(paste0(here(), "/results/para_1000_", "2021-12-03", '.csv'))
res_hal_500 <- read.csv(paste0(here(), "/results/hal_500_", "2021-12-04", '.csv'))
res_hal_1000 <- read.csv(paste0(here(), "/results/hal_1000_", "2021-12-05", '.csv'))
res_sl_500 <- read.csv(paste0(here(), "/results/sl_500_", "2021-12-04", '.csv'))
res_sl_1000 <- read.csv(paste0(here(), "/results/sl_1000_", "2021-12-05", '.csv'))

# plot
library(ggplot2)
library(ggpubr)
theme_nice <- theme(panel.background = element_blank(), 
                    panel.grid = element_line(size = 0.3, linetype = 'dashed',
                                              colour = "grey"))
library(RColorBrewer)
color.pal <- brewer.pal(n = 10, name = "Paired")

# dataset for visualization
visual_data <- rbind(tidy_visual(res_para_500), tidy_visual(res_sl_500),
                     tidy_visual(res_hal_500), tidy_visual(res_para_1000),
                     tidy_visual(res_sl_1000), tidy_visual(res_hal_1000)
                     ) %>% 
  dplyr::mutate(simu_set = rep(factor(c("param","SL","HAL")) %>% 
                                 forcats::fct_relevel("param","SL","HAL"), each = 6, 2),
                sample_size = rep(c(500, 1000), each = 18)) %>% 
  dplyr::group_by(simu_set, sample_size)

# plot of MSE and bias
plot_MSE <- visual_data %>% 
  ggplot() + 
  geom_point(aes(x= relativeMSE, y =method, shape= method), size= 4) +
  scale_shape_manual(values=rep(7,7), guide="none") + 
  geom_vline(aes(xintercept = 1), lty = 2,size=0.5) +
  geom_point(aes(x= relativeBias2, y = method, color = method), 
             shape = 21, fill = "black", size= 2)+
  scale_y_discrete(limits=rev) + 
  labs(y=NULL, x="Relative MSE and Bias^2") +
  geom_segment(aes(x = 0, xend = relativeMSE,
                   y = method, yend = method, col = method)) +
  theme_nice + facet_grid(sample_size ~ simu_set, labeller = "label_both")

# plot of coverage
plot_coverage <- visual_data %>% 
  ggplot() + 
  geom_point(aes(x= coverage, y =method, shape= method), size= 4)+
  scale_shape_manual(values=rep(7,7), guide="none") + 
  geom_vline(aes(xintercept = 0.95), lty = 2,size=0.5) +
  scale_y_discrete(limits=rev) + 
  labs(y=NULL, x="Coverage of 95% CI", legend) +
  geom_segment(aes(x = coverage-1.96*sqrt(coverage*(1-coverage)/500), 
                   xend = coverage+1.96*sqrt(coverage*(1-coverage)/500),
                   y = method, yend = method, col = method),
               arrow = arrow(angle=90, ends = "both", length = unit(0.05,"inches"))) + 
  theme_nice + facet_grid(sample_size ~ simu_set, labeller = "label_both")

# plot of relative CIwidth
plot_CIwidth <- visual_data %>% 
  ggplot() + 
  geom_point(aes(x= relativeCIwidth, y =method,  shape= method), size= 4) +
  scale_shape_manual(values=rep(7,7), guide="none") + 
  geom_vline(aes(xintercept = 1), lty = 2,size=0.5) +
  scale_y_discrete(limits=rev) + 
  theme_nice + 
  labs(y=NULL, x="Relative width of 95% CI") + 
  geom_segment(aes(x = Q1CIwidth/CIwidth * relativeCIwidth, 
                   xend = Q3CIwidth/CIwidth * relativeCIwidth,
                   y = method, yend = method, col = method),
               arrow = arrow(angle=90, ends = "both", length = unit(0.05,"inches"))) +
  theme(legend.key.width=unit(1,"cm")) +
  facet_grid(sample_size ~ simu_set, labeller = "label_both")



