# This script generates Sup Table S7 in 
# Liu et al. (2022).

# The table output shows signatures ordered by decreasing prevalence
# in the four main data sets, along with the number of runs
# of mSigHdp_ds_3k or mSigHdp_ds and SigProfilerExtractor
# in which the signature was not discovered.

if (!file.exists("common_code/set_top_dir.R")) {
  stop("Please run from top level directory")
}

require(dplyr)
require(ICAMS) # remotes::install_github("steverozen/ICAMS", ref = "v3.0.6-branch")
require(mSigTools)
require(xlsx)

source("common_code/outpath.R")

# 2. Specify global variables -------------------------------------------------

indel_or_SBS_es <- c("indel", "SBS")
folder_names <- 
  paste0(rep(indel_or_SBS_es, each = 2), 
         "_set", 
         rep(c(1L, 2L), 2))
ID_or_SBS_es <- c("ID", "SBS")
names(ID_or_SBS_es) <- indel_or_SBS_es
tool_names_general <- c("mSigHdp", "SigProfilerExtractor")
source("common_code/all.seeds.R")
seeds_in_use <- all.seeds()


# 3. Rank ground-truth signature by rareness ----------------------------------

# Read table which records the proportion of tumors with each signature
ls_dfs <- list()
# Read sig_rareness for indel_set1 and SBS_set1
sig_rareness_1 <- xlsx::read.xlsx(
  "common_code/sup_tables_s1_s2/sup_tables_s1_s2.xlsx",
  sheetIndex = 1) %>% 
  select(NA., All_type) %>%
  mutate(sig_names = NA., .keep = "unused", .before = "All_type") %>%
  mutate(sigs_prev_prop = All_type, .keep = "unused")
ls_dfs$indel_set1 <- sig_rareness_1 %>%
  filter(grepl("ID", sig_names)) %>%
  mutate(dataset_name = "indel_set1", .before = "sig_names") %>%
  arrange(desc(sigs_prev_prop))
ls_dfs$SBS_set1 <- sig_rareness_1 %>%
  filter(grepl("SBS", sig_names)) %>%
  mutate(dataset_name = "SBS_set1", .before = "sig_names") %>%
  arrange(desc(sigs_prev_prop))
# Read sig_rareness for indel_set2 and SBS_set2
sig_rareness_2 <- xlsx::read.xlsx(
  "common_code/sup_tables_s1_s2/sup_tables_s1_s2.xlsx",
  sheetIndex = 2) %>% 
  select(NA., All_type) %>%
  mutate(sig_names = NA., .keep = "unused", .before = "All_type") %>%
  mutate(sigs_prev_prop = All_type, .keep = "unused")
ls_dfs$indel_set2 <- sig_rareness_2 %>%
  filter(grepl("ID", sig_names)) %>%
  mutate(dataset_name = "indel_set2", .before = "sig_names") %>%
  arrange(desc(sigs_prev_prop))
ls_dfs$SBS_set2 <- sig_rareness_2 %>%
  filter(grepl("SBS", sig_names)) %>%
  mutate(dataset_name = "SBS_set2", .before = "sig_names") %>%
  arrange(desc(sigs_prev_prop))




# 4. Check whether each tool can extract each signature for all runs ----------
# 4a. Initialize two columns named mSigHdp and SigProfilerExtractor ------
for(fn in folder_names){
  for (tn in tool_names_general) {
    mat <- matrix(data = -1, nrow = nrow(ls_dfs[[fn]]), ncol = 1)
    colnames(mat) <- tn
    ls_dfs[[fn]] <- data.frame(ls_dfs[[fn]], mat)
  }
}

# 4b. Load Rdata object for Supp Table S4 ------
# We will fetch number of runs missed each signature from this table.
load(outpath("supplementary_table_s4.Rdata"))
supplementary_table_s4 <- supplementary_table_s4 %>%
  filter(Data_set %in% folder_names) %>%
  filter(Noise_level == "Realistic") %>%
  filter(Approach %in% c(tool_names_general, "mSigHdp_ds_3k"))


for (fn in folder_names) {
  message("fn == ", fn)
  
  # 4c. For SBS data sets, we need to use mSigHdp_ds_3k instead of mSigHdp ------
  flag_SBS <- FALSE
  if (grepl("SBS", fn) == TRUE) flag_SBS <- TRUE
  if (flag_SBS) {
    tool_names <- c("mSigHdp_ds_3k", "SigProfilerExtractor")
  } else {
    tool_names <- tool_names_general
  }
  
  for (tn in tool_names) {
    message("tn == ", tn)
    
    flag_mSigHdp_ds_3k <- FALSE
    if (tn == "mSigHdp_ds_3k") flag_mSigHdp_ds_3k <- TRUE
    for (sig_name in ls_dfs[[fn]]$sig_names) {
      # Count the number of runs which failed to extract the ground-truth sig
      # Expect to be 0~5.
      message("sig_name == ", sig_name)
      num_runs_w_fn <- 0
      for (seed_in_use in seeds_in_use) {
        message("seed_in_use ==", seed_in_use)
        stopifnot(NA %in% ls_dfs[[fn]]$sig_names == FALSE)
        message("ls_dfs[[fn]]$sig_names == ", 
                paste(ls_dfs[[fn]]$sig_names, collapse = " "))

        # Fetch 5 runs on a specific data set with a specific approach
        curr_run <- supplementary_table_s4 %>%
          filter(Data_set == fn) %>%
          filter(Approach == tn) %>%
          filter(Run == paste0("seed.", seed_in_use))
        
        FN.sigs <- unlist(curr_run$FN.sigs)
        message("FN.sigs == ", paste(FN.sigs, collapse = " "))
        
        
        if ((sig_name %in% FN.sigs) == TRUE) {
          num_runs_w_fn <- num_runs_w_fn + 1
        }
      }
      message("num_runs_w_fn ==", num_runs_w_fn) 
      row_index <- which(ls_dfs[[fn]]$sig_names == sig_name)
      if (flag_mSigHdp_ds_3k == FALSE) {
        ls_dfs[[fn]][row_index, tn] <- num_runs_w_fn
      } else {
        ls_dfs[[fn]][row_index, "mSigHdp"] <- num_runs_w_fn
      }
    }
  }
}  


# 5. Output tables ------------------------------------------------------------
# Change the order
ls_dfs <- ls_dfs[c("indel_set1", "indel_set2", "SBS_set1", "SBS_set2")]
combined_df <- do.call(rbind, ls_dfs)

utils::write.csv(combined_df, 
                 file = outpath("sup_table_s7.csv"),
                 row.names = F)



