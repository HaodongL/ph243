library(here)
source(paste0(here(), "/0_config.R"))
source(paste0(here(), "/1_hal_undersmooth.R"))
source(paste0(here(), "/2_estimation_function.R"))

##############################################################

# Part.1 Simulation

##############################################################

# ----------------------------------------
# simulation with parametric models
# ----------------------------------------
temp <- run_simu(f_simu = simu_para, 
                 psi_true = 0.1153, 
                 n_sample = 500, 
                 N_round = 2)

# ----------------------------------------
# real data based simulation with undersmoothed HAL
# ----------------------------------------

# ----------------------------------------
# real data based simulation with SL
# ----------------------------------------






###############################

# Part.2 Evaluation

###############################