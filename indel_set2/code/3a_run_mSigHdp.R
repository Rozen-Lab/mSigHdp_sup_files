# Please run this script from the top directory
if (!file.exists("common_code/set_top_dir.R")) {
  stop("Please run from top level directory")
}

message(Sys.time(), " running ", paste(commandArgs(), collapse = " "))

require(ICAMS)   # 3.0.6
require(hdpx)    # 1.0.5
require(mSigHdp) # 2.1.0 or 2.1.2 (functionally the same)


# Import optional trailing args ------------------------------------------------
curr_args <- commandArgs(trailing = T)
message("args: ", as.character(curr_args))
if (length(curr_args)==0) {
  # In this case, use the DEFAULT seed numbers
  args_flag <- FALSE
} else{
  args_flag <- TRUE
  seeds_in_use <- as.integer(curr_args)
}


# Specify global variables ----------------------------------------------------

home_for_data <- "./indel_set2/input"
home_for_run <- "./indel_set2/raw_results"

# Guessed signatures.
# We assume mSigHdp does not know the ground-truth K (13),
# then we should specify start_K as 26.
start_K <- 26

# Names of data sets to run
dataset_names <- "Realistic"

# If seeds_in_use is not specified, 
# specify 5 seeds used in software running
if (args_flag == FALSE) {
  source("common_code/all.seeds.R")
  seeds_in_use <- all.seeds()
}


# Run mSigHdp -----------------------------------------------------------------

for (dataset_name in dataset_names) {
  for (seed_in_use in seeds_in_use) {
    
    # dot case ".results" is used for compatibility with SynSigEval.
    out_dir <- paste0(home_for_run, "/mSigHdp.results/",
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
          gamma.beta      = 50,
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
    # This ensures the measurement of peak-memory usage
    # for the next job to be accurate.
    rm(multi.chains.etc)
  }
}
