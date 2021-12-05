library(here)
source(paste0(here(), "/0_config.R"))
source(paste0(here(), "/1_hal_undersmooth.R"))
source(paste0(here(), "/2_estimation_function.R"))

# parallel set up
registerDoFuture()
nCoresPerNode <- floor(as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE"))/2)
nodeNames <-strsplit(Sys.getenv("SLURM_NODELIST"), ",")[[1]]
workers <- rep(nodeNames, each=nCoresPerNode)
cl = makeCluster(workers, type = "SOCK")

plan(cluster, workers = cl)


# registerDoFuture()
# plan(multisession, workers=floor(2))

# read in data
df <- read.csv(paste0('~/Repo/ph243/results/', "washb_clean", '.csv'))

# read in hal fit
res_hal <- readRDS(paste0(here(), "/results/res_hal.RDS"))

# simulate the data
N_round = 500
n_sample = 500

df_list <- list()
for (i in 1:N_round){
  df_list[[i]] <- simu_nonp(n = n_sample, 
                            df = df, 
                            g_fit = res_hal$g_fit_hal, 
                            Q_fit = res_hal$Q_fit_hal, 
                            rv = res_hal$rv_hal , 
                            model_type = "hal")
}

# run simu
res <- run_simu(psi_true = res_hal$true_ate_hal, 
                n_sample = n_sample, 
                N_round = N_round,
                model_type = "hal",
                df_list = df_list)

output_filename <- paste0('~/Repo/ph243/results/', "hal_", n_sample, "_", Sys.Date(), '.csv')
write.csv(res, output_filename, row.names = FALSE)

stopCluster(cl)