# Please run this script from the top directory
if (basename(getwd()) != "mSigHdp_paper_sup_files_x1") {
  stop("Please run from top level directory, mSigHdp_paper_sup_files_x1")
}

# 1. Copy "Realistic" data set files under SBS_set1 ---------------------------
old_SBS_home <- "SBS_set1/input/Realistic"
SBS_home <- "SBS_set1_down_samp/input"
dir.create(SBS_home, showWarnings = FALSE, recursive = TRUE)
file.copy(old_SBS_home, SBS_home, recursive = T)
