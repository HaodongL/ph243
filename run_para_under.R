
source(paste0("~/Repo/ph243/0_config.R"))
source(paste0("~/Repo/ph243/1_hal_undersmooth.R"))
source(paste0("~/Repo/ph243/7_misQ_under.R"))

# devtools::load_all("/Users/haodongli/Repo/TMLEbootstrap/R")

# devtools::install_local("/Users/haodongli/Repo/TMLEbootstrap")
# library(TMLEbootstrap)

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
n_sample = 500

# run simu
res <- run_simu(psi_true = 0.1153, 
                n_sample = n_sample, 
                N_round = N_round,
                model_type = "para",
                df_list = NULL)


output_filename <- paste0('~/Repo/ph243/results/', "para_under_", n_sample, "_", Sys.Date(), '.csv')
write.csv(res, output_filename, row.names = FALSE)

stopCluster(cl)