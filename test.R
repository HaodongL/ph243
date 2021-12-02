library(here)
source(paste0(here(), "/0_config.R"))
source(paste0(here(), "/1_hal_undersmooth.R"))
source(paste0(here(), "/2_estimation_function.R"))

registerDoFuture()
nCoresPerNode <- floor(as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE"))/4)
nodeNames <-strsplit(Sys.getenv("SLURM_NODELIST"), ",")[[1]]
workers <- rep(nodeNames, each=nCoresPerNode)
cl = makeCluster(workers, type = "SOCK")

plan(cluster, workers = cl)


res <- run_simu(f_simu = simu_para, 
                 psi_true = 0.1153, 
                 n_sample = 500, 
                 N_round = 4)


output_filename <- paste0('~/Repo/ph243/results/', "test_", Sys.Date(), '.csv')
write.csv(res, output_filename, row.names = FALSE)

stopCluster(cl)
