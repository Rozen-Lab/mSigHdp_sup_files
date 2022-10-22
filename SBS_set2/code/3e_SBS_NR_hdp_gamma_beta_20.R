# How to run:
# $ nice Rscript SBS_set2/code/4e_SBS_run_NR_hdp_gamma_beta_20.R <seed_number> >>& SBS_set2/raw_results/NR_hdp_gb_20.results/log &

# Please run this script from the top directory
if (!file.exists("common_code/set_top_dir.R")) {
  stop("Please run from top level directory")
}


message(Sys.time(), " running ", paste(commandArgs(), collapse = " "))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  seeds_in_use <- args
} else {
  source("common_code/all.seeds.R")
  seeds_in_use <- all.seeds()
}

message(Sys.time(), " running on seed ", seeds_in_use)

# Set global variables ---------------------------------------------------------
GLOBAL.gamma.alpha <- 1  # This will be used inside mSigHdp::SetupAndPosterior;
                         # alpha is also called the shape parameter
GLOBAL.gamma.beta  <- 20 # This will be used inside mSigHdp::SetupAndPosterior;
                         # beta is also called the rate parameter
burnin.iterations  <- 5000 * 20
CPU.cores          <- 20
num.child.process  <- 20
# Guessed number of raw clusters
start_K            <- 64

home_for_run       <- paste0("./SBS_set2/raw_results/NR_hdp_gb_", 
                             GLOBAL.gamma.beta, ".results/")
home_for_data      <- "./SBS_set2/input"

# Names of data sets
dataset_names <- c("Realistic")

# Run mSigHdp -----------------------------------------------------------------
 
source("common_code/generic_run_NR_hdp.R")
