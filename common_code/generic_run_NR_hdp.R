if (!file.exists("common_code/set_top_dir.R")) {
  stop("Please run from top level directory")
}

if (packageVersion("hdpx") != "0.1.5.0099" ||
  packageVersion("mSigHdp") != "0.0.0.9019") {
  stop(
    "To test Nicola Roberts's algorithm for combining raw mutation\n",
    "clusters, you must use specific versions of hdpx and mSigHdp.\n",
    "Please run common_code/install_NR_hdp.R to install these.")
}

library(ICAMS)
library(hdpx)
library(mSigHdp)

message(Sys.time(), " starting analyses; output in subdir of ", home_for_run)
if (!dir.exists(home_for_run)) {
  message(Sys.time(), " creating ", home_for_run)
  dir.create(home_for_run, recursive = TRUE)
} else {
  message(Sys.time(), home_for_run, " already exists")
}

# Run mSigHdp -----------------------------------------------------------------

for (dataset_name in dataset_names) {
  for (seed_in_use in seeds_in_use) {
    
    # dot case ".results" is used for compatibility with SynSigEval.
    out_dir <- paste0(home_for_run, "/",
                      dataset_name, "/seed.", seed_in_use)
    
    message(Sys.time(), " putting results in ", out_dir)
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    }
    
    seed_in_use <- as.numeric(seed_in_use)
    message(Sys.time(), " seed_in_use=", seed_in_use)

    
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
        multi.chains.etc <- mSigHdp::RunHdpParallel(
          input.catalog     = paste0(home_for_data, "/", dataset_name,
                                     "/ground.truth.syn.catalog.csv"),
          seedNumber        = seed_in_use,
          K.guess           = start_K,
          multi.types       = TRUE,
          post.burnin       = burnin.iterations,
          post.n            = 200,
          post.space        = 100,
          CPU.cores         = CPU.cores,
          num.child.process = num.child.process,
          overwrite         = TRUE,
          out.dir           = out_dir)
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
