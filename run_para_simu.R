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
print(paste0("nCoresPerNode: ", nCoresPerNode))

# registerDoFuture()
# plan(multisession, workers=floor(2))

N_round = 500
n_sample = 1000

# run simu
res <- run_simu(psi_true = 0.1153, 
                n_sample = n_sample, 
                N_round = N_round,
                model_type = "para",
                df_list = NULL)


output_filename <- paste0('~/Repo/ph243/results/', "para_", n_sample, "_", Sys.Date(), '.csv')
write.csv(res, output_filename, row.names = FALSE)

stopCluster(cl)
