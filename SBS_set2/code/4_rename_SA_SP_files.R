# Please run this script from the top directory
if (!file.exists("common_code/set_top_dir.R")) {
  stop("Please run from top level directory")
}

require(ICAMS)   # 3.0.6



# Specify global variables ----------------------------------------------------
home_for_data <- "SBS_set2/input"
home_for_run <- "SBS_set2/raw_results"

# Import data set names
dataset_names <- "Realistic"

# Specify 5 seeds used in software running
seeds_in_use <- c(145879, 200437, 310111, 528401, 1076753)



# Copy the best extraction by SignatureAnalyzer in "best.run/" ----------------
for (dataset_name in dataset_names) {
  for (seed_in_use in seeds_in_use) {
    
    sa_run_dir <- paste0(home_for_run, "/SignatureAnalyzer.results/",
                         dataset_name, "/seed.", seed_in_use)
    
    # Copy signature in ICAMS format
    sig_path <- paste0(sa_run_dir, "/best.run/sa.output.sigs.csv")
    file.copy(from = sig_path,
              to = paste0(sa_run_dir, "/extracted.signatures.csv"),
              copy.date = TRUE)
  }
}

# Copy SigPro signature files from old.sig.path to sig.path -------------------
for (dataset_name in dataset_names) {
  for (seed_in_use in seeds_in_use) {
    sp_run_dir <- paste0(home_for_run, "/SigProfilerExtractor.results/",
                         dataset_name, "/seed.", seed_in_use)
    # Note: SigProfilerExtractor originally exports extracted signature file
    # in a very deep path
    old.sig.path <- 
      paste0(sp_run_dir, 
             "/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Signatures/",
             "SBS96_De-Novo_Signatures.txt")
    # This very deep path will break some program on Windows,
    # and thus we copied these files to sig.path.
    sig.path <- paste0(sp_run_dir, "/SBS96_De-Novo_Signatures.txt")
    file.copy(from = old.sig.path, to = sig.path, copy.date = TRUE)
  }
}


# Convert SigPro-TSV-formatted signatures to ICAMS-CSV format -----------------
for (dataset_name in dataset_names) {
  for (seed_in_use in seeds_in_use) {
    sp_run_dir <- paste0(home_for_run, "/SigProfilerExtractor.results/",
                         dataset_name, "/seed.", seed_in_use)
    
    # Convert catalog to ICAMS format, using wrapper function
    sig.path <- paste0(sp_run_dir, "/SBS96_De-Novo_Signatures.txt")
    sig.catalog.sp <- utils::read.table(
      sig.path,
      sep = "\t",
      as.is = TRUE,
      header = TRUE)
    sig.catalog <- ICAMS:::MakeSBS96CatalogFromSigPro(sig.catalog.sp)
    sig.catalog <- ICAMS::as.catalog(sig.catalog,
                                     catalog.type = "counts.signature")
    ICAMS::WriteCatalog(sig.catalog,
                        paste0(sp_run_dir, "/extracted.signatures.csv"))
  }
}
