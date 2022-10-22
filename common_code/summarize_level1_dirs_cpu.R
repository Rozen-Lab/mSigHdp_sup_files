# To include SigProfilerExtractor's CPU time
# into the summary, please install pandas
# in the python environment configured in
# reticulate::use_python().


if (!file.exists("common_code/set_top_dir.R")) {
  stop("Please run from top level directory")
}

require(openxlsx)
require(reticulate)
require(tibble)

source("common_code/outpath.R")

# Need to have pandas installed in your python environment
#reticulate::py_install("pandas", ignore_installed = FALSE)
pd <- reticulate::import("pandas")



# Specify global variables ----------------------------------------------------
# Level 1 dirs are also the name of the data sets
level1_dirs <- c("indel_set1",
                 "indel_set2",
                 "SBS_set1",
                 "SBS_set2")
level1_dirs <- c(level1_dirs, "SBS_set1_down_samp", "SBS_set2_down_samp")

if (FALSE) {
level1_dirs <- c(level1_dirs, paste0(level1_dirs, "_down_samp"))
level1_dirs <- c(level1_dirs,
                 paste0("sens_SBS35_",
                        c(5L, 10L, 20L, 30L, 50L, 100L),
                        "_1066"))
}


# Function for summary directory on one seed ----------------------------------
cpu_time_from_one_seed <- 
  function(run_directory_path, SP_flag = FALSE) {
  SP_flag <- FALSE
  if (grepl("SigProfilerExtractor", run_directory_path)) SP_flag <- TRUE
  if (grepl("SP_narrower_range", run_directory_path)) SP_flag <- TRUE
  
  if(SP_flag == TRUE) {
    # Load profiling_info object from python pickle file
    profiling_info <- pd$read_pickle(
      file.path(run_directory_path, "profiling_info.pickle")
    )
    cpu_time <- profiling_info$cpu_time
  } else {
    # Load code.profile object
    load(file.path(run_directory_path, "code.profile.Rdata"))
    # CPU time in seconds
    user_cpu <- code.profile$system.time[c(1,4)] %>% sum()
    system_cpu <- code.profile$system.time[c(2,5)] %>% sum()
    cpu_time <- user_cpu + system_cpu
  }
  return(cpu_time)
}

cpu_time_level1_dirs <- function(a_folder) {
  message("summarizing a_folder=", a_folder)
  stopifnot(dir.exists(a_folder))
  dataset_name <- sub("_down_samp", "", basename(a_folder))
  message("using dataset name ", dataset_name)

  cpu_time_table <- tibble_row(Data_set         = "",
                               Approach         = "",
                               Run              = "",
                               cpu_time         = -1)
  start_here <- file.path(a_folder, "raw_results")
  
  tool_paths <- list.files(start_here, 
                           pattern = "\\.results$",
                           full.names = TRUE)
  # Remove "SP_KmSigHdp.results" and "SP_Kplus2.results"
  tool_paths <- 
    tool_paths[!(basename(tool_paths) %in% 
               c("SP_KmSigHdp.results", "SP_Kplus2.results"))]
  
  for (tool_path in tool_paths) {
    # tool_path - a full directory path, 
    # e.g "indel/raw_results/mSigHdp.results"
    stopifnot(dir.exists(tool_path))
    
    message("Checking CPU times of tool_path=", tool_path)
    tool_name <- sub(".results", "", basename(tool_path), fixed = TRUE)
    
    if (!tool_name %in% c("mSigHdp",
                          "mSigHdp_ds_3k",
                          "SigProfilerExtractor",
                          "SignatureAnalyzer",
                          "signeR")) next
    
    dataset_paths <- list.files(tool_path, full.names = TRUE)
    for(dataset_path in dataset_paths){
      noise_level <- basename(dataset_path)
      if (!noise_level =="Realistic") {
        if (dir.exists(dataset_path)) {
          message("\n\n**Skipping directory ", dataset_path, "**\n\n")
        }
        next
      }
      message("Checking CPU times of dataset_path=", dataset_path)
      seeds_paths <- list.files(dataset_path, 
                                pattern = "seed\\.\\d+$", 
                                full.names = TRUE)
      
      for(seed_path in seeds_paths) {
        # seed_path e.g. "indel/raw_results/mSigHdp.results/Moderate/seed.528401"
        message("Checking CPU time for seed_path=", seed_path)
        cpu_time <- cpu_time_from_one_seed(seed_path)
        a_row <- tibble_row(Data_set         = dataset_name,
                            Approach         = tool_name,
                            Run              = basename(seed_path),
                            cpu_time         = cpu_time)

        cpu_time_table <- rbind(cpu_time_table, a_row)
      } # for(seed_path in seeds_paths)
    } # for (dataset_path in dataset_names)
  } # for (tool_path in tool_paths)
  
  # Remove the first dummy row
  cpu_time_table <- cpu_time_table[-1, ]
  # readr::write_csv(
  #   data.table::as.data.table(cpu_time_table), 
  #   file.path(a_folder, "level1_cpu_time.csv"))
  # save(cpu_time_table, file=file.path(a_folder, "level1_cpu_time.Rdata"))
  invisible(return(cpu_time_table))
} # function cpu_time_something


# Wrapper function ------------------------------------------------------------
cpu_time_all_level1_dirs <- function(level1_dirs = level1_dirs) {

  all_out_list <- lapply(level1_dirs, cpu_time_level1_dirs)
  
  all_cpu_time <- do.call(rbind, all_out_list)
  
  utils::write.csv(all_cpu_time, 
                   outpath("supplementary_table_s5.csv"), 
                   row.names = F,
                   quote = F)
  openxlsx::write.xlsx(all_cpu_time, 
                       outpath("supplementary_table_s5.xlsx"))
  save(all_cpu_time, file = outpath("supplementary_table_s5.Rdata"))
  invisible(all_cpu_time)
}

# Run wrapper function --------------------------------------------------------
all_cpu_time <- invisible(cpu_time_all_level1_dirs(level1_dirs))
