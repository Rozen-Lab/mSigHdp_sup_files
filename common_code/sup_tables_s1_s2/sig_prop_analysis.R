# Please run this script from the top directory
if (!file.exists("common_code/set_top_dir.R")) {
  stop("Please run from top level directory")
}

library(mSigTools)
library(PCAWG7) # remotes::install_github(repo = "steverozen/PCAWG7", ref = "v0.1.3-branch")
library(dplyr)
library(gtools)
library(xlsx)

sbs_exp1 <-
  mSigTools::read_exposure(file = "SBS_set1/input/Realistic/ground.truth.syn.exposures.csv")
sbs_exp2 <-
  mSigTools::read_exposure(file = "SBS_set2/input/Realistic/ground.truth.syn.exposures.csv")

indel_exp1 <-
  mSigTools::read_exposure(file = "indel_set1/input/Realistic/ground.truth.syn.exposures.csv")
indel_exp2 <-
  mSigTools::read_exposure(file = "indel_set2/input/Realistic/ground.truth.syn.exposures.csv")

get_sig_prop <- function(exp1) {
  get_prop <- function(exposure, identifier) {
    prop <-
      apply(X = exposure, MARGIN = 1, FUN = function(x) {
        length(x[x > 0]) / length(x)
      })
    df <- as.data.frame(t(prop))
    rownames(df) <- identifier
    return(df)
  }
  all_prop1 <- get_prop(exposure = exp1, identifier = "All_type")

  exp_by_type1 <-
    PCAWG7::SplitPCAWGMatrixByTumorType(exp1)

  cancer_types1 <-
    gsub(pattern = "SP.Syn.", replacement = "", x = names(exp_by_type1))

  cancer_types <- cancer_types1

  prop_list <- list(all_prop1)

  for (cancer_type in cancer_types) {
    exp1_samples <-
      grep(pattern = cancer_type, x = colnames(exp1), value = TRUE)
    exp_type1 <- exp1[, exp1_samples]
    exp_type1 <- exp_type1[rowSums(exp_type1) > 0, , drop = FALSE]
    prop1 <- get_prop(exposure = exp_type1, identifier = cancer_type)
    prop_list <- c(prop_list, list(prop1))
  }

  prop_df <- do.call(dplyr::bind_rows, prop_list)
  prop_df2 <- as.data.frame(t(prop_df))
  prop_df3 <- prop_df2[gtools::mixedsort(x = rownames(prop_df2)), , drop = FALSE]
  return(prop_df3)
}

sbs_prop_set1 <- get_sig_prop(exp1 = sbs_exp1)
indel_prop_set1 <- get_sig_prop(exp1 = indel_exp1)

sbs_prop_set2 <- get_sig_prop(exp1 = sbs_exp2)
indel_prop_set2 <- get_sig_prop(exp1 = indel_exp2)

sig_prop_set1 <- rbind(sbs_prop_set1, indel_prop_set1)
sig_prop_set2 <- rbind(sbs_prop_set2, indel_prop_set2)
xlsx::write.xlsx(sig_prop_set1,
  file = "common_code/sup_tables_s1_s2/sup_tables_s1_s2.xlsx",
  showNA = FALSE, sheetName = "set1"
)

xlsx::write.xlsx(sig_prop_set2,
  file = "common_code/sup_tables_s1_s2/sup_tables_s1_s2.xlsx",
  showNA = FALSE, sheetName = "set2",
  append = TRUE
)
