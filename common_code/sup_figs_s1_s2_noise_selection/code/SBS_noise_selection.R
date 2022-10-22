# Please run this script from the top directory
if (!file.exists("common_code/set_top_dir.R")) {
  stop("Please run from top level directory")
}

source("./common_code/data_gen_utils.R")
source("common_code/outpath.R")

library(dplyr)
library(ggpubr)
library(gridExtra)
library(cosmicsig)
library(ICAMS) # remotes::install_github("steverozen/ICAMS", ref = "v3.0.6-branch")
library(PCAWG7) # remotes::install_github(repo = "steverozen/PCAWG7", ref = "v0.1.3-branch")
library(mSigAct) # remotes::install_github(repo = "steverozen/mSigAct", ref = "v2.2.0-branch")
library(SynSigGen) # remotes::install_github(repo = "steverozen/SynSigGen", ref = "1.1.1-branch")

# Get the real exposure for tumors from selected cancer types
real_exposure_sbs_file <-
  "./common_code/sup_figs_s1_s2_noise_selection/data/SBS_real_exposure.csv"
real_exposure_sbs <- mSigAct::ReadExposure(file = real_exposure_sbs_file)
sigs_sbs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS96

# Get the real tumor spectra from the selected cancer types
real_spectra_sbs <-
  PCAWG7::spectra$PCAWG$SBS96[, colnames(real_exposure_sbs)]

real_distance_sbs <- get_distance(
  spectra = real_spectra_sbs,
  exposure = real_exposure_sbs,
  sigs = sigs_sbs,
  group = "real"
)

# Add noise to noiseless synthetic data with different negative-binomial size
# parameter
dir_noiseless_sbs <- "./SBS_set2/input/Noiseless/"
noiseless_data_sbs <- get_syn_data_info(dir = dir_noiseless_sbs)
seed <- 658220

nb_sizes_sbs <- c(100, 50, 40, 30, 20, 10)
noise_data_sbs <- generate_noisy_data(
  seed = seed,
  exposure = noiseless_data_sbs$exposure,
  sigs = noiseless_data_sbs$sigs,
  nb_sizes = nb_sizes_sbs
)
syn_distance_sbs <-
  get_multiple_syn_distances(list_of_syn_data = noise_data_sbs)

all_distance_sbs <-
  do.call(dplyr::bind_rows, c(list(real_distance_sbs), syn_distance_sbs))

plot_objects_sbs <-
  multiple_boxplots(
    distance_df = all_distance_sbs,
    data_type = "SBS", ylim = c(0, 0.2)
  )

ggplot_to_pdf(
  plot_objects = plot_objects_sbs,
  file = outpath("sup_fig_s1_SBS_noise_selection.pdf"),
  nrow = 2, ncol = 1,
  width = 8.2677, height = 11.6929, units = "in"
)

