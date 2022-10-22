# Please run this script from the top directory
if (basename(getwd()) != "Liu_et_al_Sup_Files") {
  stop("Please run from top level directory, Liu_et_al_Sup_Files")
}

# 1. Specify global options ---------------------------------------------------
options(stringsAsFactors = F)



# 2. Specify global variables -------------------------------------------------

# Specify source and destination folder for file copy.
old_home_for_run <- "./indel/raw_results/"
new_home_for_run <- "./other_analyses/indel_beta_selection/raw_results"

# Specify name of ground-truth spectra data set
dataset_name <- "Realistic"

# Specify names of computational approach
tool_name <- "SigProfilerExtractor"

# Specify seeds used in analysis.
# Specify 5 seeds used in software running
seeds_in_use <- c(145879, 200437, 310111, 528401, 1076753)

# Specify file name to be copied.
file_name <- c("extracted.signatures.csv")


# 3. Create sub-directories under destination folder --------------------------

for (seed_in_use in seeds_in_use) {
  dir.create(paste0(new_home_for_run, "/", 
                    tool_name, ".results/",
                    dataset_name, "/",
                    "seed.", seed_in_use, "/"),
             recursive = T)
}



# 4. Copy SigProfilerExtractor files ------------------------------------------

# Copies signatures extracted by SigProfilerExtractor 
# from Realistic data sets.
for (seed_in_use in seeds_in_use) {
  source_file <- paste0(old_home_for_run, "/", 
                        tool_name, ".results/",
                        dataset_name, "/",
                        "seed.", seed_in_use, "/", file_name)
  dest_file <- paste0(new_home_for_run, "/", 
                      tool_name, ".results/",
                      dataset_name, "/",
                      "seed.", seed_in_use, "/", file_name)
  file.copy(source_file, dest_file, copy.date = T)
}
