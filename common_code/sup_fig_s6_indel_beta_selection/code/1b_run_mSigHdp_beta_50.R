# Please run this script from the top directory
if (basename(getwd()) != "Liu_et_al_Sup_Files") {
  stop("Please run from top level directory, Liu_et_al_Sup_Files")
}

# Install and load required packages ------------------------------------------

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
# "@*release" installs the latest release
if (!("hdpx" %in% rownames(installed.packages())) ||
    packageVersion("hdpx") < "1.0.1") {
  remotes::install_github("steverozen/hdpx@*release")
}
if (!("mSigHdp" %in% rownames(installed.packages())) ||
    packageVersion("mSigHdp") < "2.0.0") {
  remotes::install_github(repo = "steverozen/mSigHdp@*release")
}

# ICAMS is installed when installing mSigHdp
require(ICAMS)
require(hdpx)
require(mSigHdp)



# Specify global variables ----------------------------------------------------

home_for_data <- "./indel/input"
home_for_run <- "./other_analyses/indel_beta_selection/raw_results"

# Guessed signatures.
# We assume mSigHdp does not know the ground-truth K (11),
# then we should specify start_K as 22.
start_K <- 22

# Value of gamma.beta to use
gamma.beta <- 50

# Names of data set
dataset_names <- "Realistic"

# Specify 5 seeds used in software running
seeds_in_use <- c(145879, 200437, 310111, 528401, 1076753)



# Run mSigHdp -----------------------------------------------------------------

for (dataset_name in dataset_names) {
  for (seed_in_use in seeds_in_use) {
    
    # dot case ".results" is used for compatibility with SynSigEval.
    out_dir <- paste0(home_for_run,"/mSigHdp.beta", gamma.beta,".results/",
                      dataset_name, "/seed.", seed_in_use)
    
    # Skip if all finished jobs to save time if a users needs to re-run.
    if (file.exists(paste0(out_dir, "/code.profile.Rdata"))) next
    
    message("\n\n===========================================\n\n")
    message(paste0("Begin running mSigHdp on data set ",
                   dataset_name, " using seed ", seed_in_use, "...\n"))
    message("\n\n===========================================\n\n")
    
    # Instantiate a list to store profiling data.
    code.profile <- list()
    
    # Garbage collection before code profiling 
    # to make code profiling more accurate
    #
    # Specify reset = TRUE to reset max memory consumption.
    # This prevents the influence by previous runs. 
    gc(reset = TRUE)
    
    # system.time() is a length-5 vector which collects the following 
    # time usage:
    #
    # 1. user.self: User CPU time of the main process
    # 2. system.self: System CPU time of the main process
    # 3. elapsed: Elapsed time (a.k.a. Wall clock time) - 
    #    the actual time a user experienced to wait from the job's 
    #    beginning to the end. 
    # 4. user.child: Sum of user CPU times of all the child processes
    # 5. system.child: Sum of system CPU times of all the child processes
    #
    # When calling "print(code_profile$system.time)" only three numbers
    # are returned:
    #
    #   user: Sum of 1 and 4; 
    #   system: Sum of 2 and 5;
    #   elapsed: 3.
    code.profile[["system.time"]] <- system.time(
      {
        multi.chains.etc <- mSigHdp::RunHdpxParallel(
          input.catalog = paste0(home_for_data, "/", dataset_name,
                                 "/ground.truth.syn.catalog.csv"),
          seedNumber = seed_in_use,
          out.dir = out_dir,
          K.guess = start_K,
          multi.types = TRUE,
          # Number of burn-in iterations in each posterior chain.
          # mSigHdp::RunHdpxParallel() uses 5000 ("burnin") * 2 
          # ("burnin.multiplier") by default.
          #
          # From prior testing we found that the likelihood 
          # has not converged when using burnin.multiplier = 2.
          #
          # Therefore, we need to extend the burn-in iterations 3X,
          # Causing the total number of iterations +67% (30,000 -> 50,000)
          burnin     = 5000,
          burnin.multiplier = 6,
          post.n          = 200,
          post.space      = 100,
          num.child.process = 20,
          CPU.cores = 20,
          high.confidence.prop = 0.9,
          gamma.alpha     = 1,
          gamma.beta      = gamma.beta,
          overwrite       = T)
        save(multi.chains.etc, file = paste0(out_dir, "/multi.chains.etc.Rdata"))
      },
      gcFirst = FALSE
    )
    
    # The return value of gc()
    # records info for peak memory usage.
    code.profile[["gc"]] <- gc(reset = TRUE)
    
    # Save code profiling data.
    save(code.profile, file = paste0(out_dir, "/code.profile.Rdata"))
    
    # Delete mSigHdp return-value object.
    # This makes the measurement of peak-memory usage
    # for the next job accurate.
    rm(multi.chains.etc)
  }
}
