# Small calculations cited in the text
# Results are usually just printed to the console.

library(data.table)
library(tidyr)
library(MASS)
library(dplyr)
library(mSigHdp)

source("common_code/outpath.R")


# Ratio of number of mutations in downsample(SBS_set2) to number of -----------
# mutations in downsample(SBS_set1)

sset1 <- ICAMS::ReadCatalog(
  "SBS_set1/input/Realistic/ground.truth.syn.catalog.csv")
sset2 <- ICAMS::ReadCatalog(
  "SBS_set2/input/Realistic/ground.truth.syn.catalog.csv")

count_sset1 <- sum(sset1)
count_sset2 <- sum(sset2)
(count_sset1 + count_sset2) / 1e6

count_sset2 / count_sset1
sum(mSigHdp::downsample_spectra(sset2, downsample_threshold = 3000)[[1]]) /
  sum(mSigHdp::downsample_spectra(sset1, downsample_threshold = 3000)[[1]])


# Statistics (robust regression) for Supplementary Table S7 -------------------

s7 <- fread(outpath("sup_table_s7.csv"))
s7x <- tidyr::pivot_longer(s7, 
                           cols = c("mSigHdp", "SigProfilerExtractor"),
                           names_to = "Approach",
                           values_to = "num_seeds_with_miss")
s7x.indel <- dplyr::filter(s7x, grepl("indel", dataset_name))
s7x.sbs <- dplyr::filter(s7x, grepl("SBS", dataset_name))
sup_table_s7_sbs_for_rlm <- 
  dplyr::mutate(s7x.sbs, num_seeds_with_hit = 5 - num_seeds_with_miss)
data.table::fwrite(sup_table_s7_sbs_for_rlm, 
                   outpath("sup_table_s7_sbs_for_rlm.csv"))

sbsr <- MASS::rlm(num_seeds_with_hit ~ sigs_prev_prop + Approach + 
                    dataset_name, data = sup_table_s7_sbs_for_rlm)
rm   <- summary(sbsr)
rmc  <- rm$coefficients
rmp  <- 2*pt(-abs(rmc[ , 3]), df = nrow(sup_table_s7_sbs_for_rlm) - 4)
coef <- data.frame(cbind(rmc, p = rmp))
sup_table_s7_coef <- cbind(Variable = rownames(coef), coef)
openxlsx::write.xlsx(sup_table_s7_coef, 
                     file = outpath("sup_table_s7_coef.xlsx"))


# Supplementary Table S7, double checking numbers of false negatives for ------
# SBS_set1
load(outpath("supplementary_table_s4.Rdata"))

dplyr::filter(supplementary_table_s4,
              Approach %in% c("mSigHdp_ds_3k", "SigProfilerExtractor")) %>%
  dplyr::filter(grepl("SBS_set", Data_set)) %>%
  dplyr::filter(Noise_level == "Realistic") %>%
  dplyr::group_by(Approach) -> sup.table.s7
tidyr::unnest_longer(sup.table.s7, FN.sigs) -> st7
dplyr::filter(st7, FN.sigs == "SBS29")
dplyr::filter(st7, FN.sigs == "SBS30")


# Numbers for the section on indel results ------------------------------------


dplyr::filter(supplementary_table_s4,
                Data_set %in% c("indel_set1", "indel_set2") &
                Noise_level == "Realistic") %>%
  dplyr::filter(Approach %in% 
                  c("mSigHdp", 
                    "NR_hdp_gb_1", 
                    "NR_hdp_gb_50", 
                    "SigProfilerExtractor")) %>%
  dplyr::group_by(Approach) %>%
  dplyr::mutate(mean.cm = mean(Composite), 
                med.fp = median(FP), 
                avg.fp = mean(FP)) -> examine.indels
  View(examine.indels)
  dplyr::distinct(dplyr::transmute(examine.indels, Approach = Approach, mean.cm = mean.cm))
  

# Main text Page 11, means for SBS data across --------------------------------
# both SBS_set1 and SBS_set2

dplyr::filter(supplementary_table_s4, grepl("SBS_set", Data_set)) %>%
  dplyr::filter(Noise_level == "Realistic") %>%
  dplyr::group_by(Approach) %>%
  dplyr::transmute(Data_set = Data_set,
                   med.cm = median(Composite), 
                   avg.cm = mean(Composite),
                   avg.ppv = mean(PPV),
                   avg.tpr = mean(TPR)) -> p11.step1
dplyr::filter(p11.step1, 
              Approach %in% c("mSigHdp", "mSigHdp_ds_3k",
                              "SigProfilerExtractor",
                              "NR_hdp_gb_1", "NR_hdp_gb_20",
                              "SignatureAnalyzer", "signeR")) -> p11.step2

dplyr::mutate(p11.step2, Data_set = NULL) %>% dplyr::distinct()


# Supplementary Figure S7, comparison of average Composite Measure over -------
# SBS_set1 and SBS_set for different downsampling thresholds.

dplyr::filter(supplementary_table_s4,
              grepl("^mSigHdp", Approach)) %>%
  dplyr::filter(grepl("SBS_set", Data_set)) %>%
  dplyr::filter(Noise_level == "Realistic") %>%
  dplyr::group_by(Approach) %>%
  dplyr::transmute(Data_set = Data_set,
                   med.cm = median(Composite), 
                   avg.cm = mean(Composite),
                   avg.ppv = mean(PPV),
                   avg.tpr = mean(TPR)) -> summary.sup.fig.s7.step1

dplyr::mutate(summary.sup.fig.s7.step1, Data_set = NULL) %>%
  dplyr::distinct() -> summary.sup.fig.s7


# Find best values of K from mSigHdp for sup table s6

dplyr::filter(supplementary_table_s4, 
Data_set %in% c("SBS_set1", "indel_set1")) %>%
  dplyr::filter(Noise_level == "Realistic" &
                  Approach %in% c("mSigHdp", "mSigHdp_ds_3k")) -> foo
dplyr::group_by(foo, Approach, Data_set) %>%
  dplyr::transmute(Approach = Approach,
                   Data_set = Data_set,
                  med.N.sig = median(N_Sigs),
                   avg.N.sig = mean(N_Sigs),
                  N_Sigs     = N_Sigs,
                  Composite = Composite,
                  best.comp  = max(Composite)) -> foo2

indel_ms_K <- dplyr::filter(
  foo2, Data_set == "indel_set1" & Approach == "mSigHdp")
dplyr::filter(indel_ms_K, Composite == best.comp)


SBS_ms_K   <- 
  dplyr::filter(foo2, Data_set == "SBS_set1" & Approach == "mSigHdp_ds_3k")
dplyr::filter(SBS_ms_K, Composite == best.comp)
rm(foo, foo2)

