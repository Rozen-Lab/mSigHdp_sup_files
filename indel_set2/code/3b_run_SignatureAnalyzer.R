# Please run this script from the top directory
if (!file.exists("common_code/set_top_dir.R")) {
  stop("Please run from top level directory")
}

require(ICAMS)      # 3.0.6 (currently at github steverozen/ICAMS)
require(SynSigRun)  # 1.0.0 remotes::install_github("WuyangFF95/SynSigRun", ref = "1.0.0-branch")


# Download SynSigRun to fetch SignatureAnalyzer source code -------------------

# Do not insert trailing "/" to prevent failing unzip
SA_home <- "./common_code"

SA_path <- paste0(SA_home, "/SignatureAnalyzer.052418")

if (!dir.exists(SA_path)) {
  dir.create(SA_home, showWarnings = F, recursive = T)
  if (!file.exists(paste0(SA_path, ".zip"))) {
    download.file(
      "https://github.com/WuyangFF95/SynSigRun/raw/master/data-raw/SignatureAnalyzer.052418.zip",
      destfile = paste0(SA_home,"/SignatureAnalyzer.052418.zip")
    )
  }
  unzip(paste0(SA_home, "/SignatureAnalyzer.052418.zip"), exdir = SA_home)
  unlink(paste0(SA_home, "/SignatureAnalyzer.052418.zip"))
}



# Import optional trailing args ------------------------------------------------
curr_args <- commandArgs(trailing = T)
message("args: ", as.character(curr_args))
if (length(curr_args)==0) {
  # In this case, use the DEFAULT seed numbers
  args_flag <- FALSE
} else {
  args_flag <- TRUE
  seeds_in_use <- as.integer(curr_args)
}



# Specify global variables ----------------------------------------------------

home_for_data <- "./indel_set2/input"
home_for_run <- "./indel_set2/raw_results"

# Guessed signatures.
# We assume SignatureAnalyzer does not know the ground-truth K (13),
# then we should specify max_K as 26.
max_K <- 26

# Names of data sets
dataset_names <- c("Realistic")

# If seeds_in_use is not specified,
# specify 5 seeds used in software running
if (args_flag == FALSE) {
  source("common_code/all.seeds.R")
  seeds_in_use <- all.seeds()
}


# Run SignatureAnalyzer -------------------------------------------------------

for (dataset_name in dataset_names) {
  for (seed_in_use in seeds_in_use) {
    
    # dot case ".results" is used for compatibility with SynSigEval.
    out_dir <- paste0(home_for_run,"/SignatureAnalyzer.results/",
                      dataset_name,"/seed.",seed_in_use)
    
    # Skip if all finished jobs to save time if a users needs to re-run.
    if(file.exists(paste0(out_dir,"/code.profile.Rdata"))) next
    
    message("\n\n===========================================\n\n")
    message(paste0("Begin running SignatureAnalyzer on data set ",
                   dataset_name," using seed ",seed_in_use,"...\n"))
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
        # SignatureAnalyzer needs to run 20 parallel runs,
        # and pick the best run as the final result.
        SynSigRun::SAMultiRunOneCatalog(
          num.runs = 20,
          signatureanalyzer.code.dir = SA_path,
          input.catalog = paste0(home_for_data,"/",dataset_name,
                                 "/ground.truth.syn.catalog.csv"),
          out.dir = out_dir,
          maxK = max_K,
          tol = 1e-7,
          test.only = FALSE,
          delete.tmp.files = TRUE,
          overwrite = FALSE,
          mc.cores = 1,
          verbose = FALSE,
          seed = seed_in_use)
        
        SynSigRun:::CopyBestSignatureAnalyzerResult(
          sa.results.dir = out_dir,
          verbose = TRUE,
          overwrite = FALSE)
      },
      gcFirst = FALSE
    )
    # Garbage collection.
    # The return value of gc()
    # records info for peak memory usage.
    code.profile[["gc"]] <- gc(reset = TRUE)
    
    # Save code profiling data.
    save(code.profile, file = paste0(out_dir,"/code.profile.Rdata"))
  }
}
