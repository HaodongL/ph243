library(here)
source(paste0(here(), "/0_config.R"))
source(paste0(here(), "/1_hal_undersmooth.R"))
source(paste0(here(), "/2_estimation_function.R"))

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

# read in sl fit
res_sl <- readRDS(paste0(here(), "/results/res_sl.RDS"))

# simulate the data
N_round = 500
n_sample = 500

df_list <- list()
for (i in 1:N_round){
  df_list[[i]] <- simu_nonp(n = n_sample, 
                            df = df, 
                            g_fit = res_sl$g_fit_sl, 
                            Q_fit = res_sl$Q_fit_sl, 
                            rv = res_sl$rv_sl , 
                            model_type = "sl")
}

# run simu
res <- run_simu(psi_true = res_sl$true_ate_sl, 
                n_sample = n_sample, 
                N_round = N_round,
                model_type = "sl",
                df_list = df_list)

output_filename <- paste0('~/Repo/ph243/results/', "sl_", n_sample, "_", Sys.Date(), '.csv')
write.csv(res, output_filename, row.names = FALSE)

stopCluster(cl)