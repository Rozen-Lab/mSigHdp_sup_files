# Please run this script from the top directory
if (!file.exists("common_code/set_top_dir.R")) {
  stop("Please run from top level directory")
}

source("./common_code/data_gen_utils_set2.R")

library(cosmicsig) 
library(PCAWG7) # remotes::install_github(repo = "steverozen/PCAWG7", ref = "v0.1.3-branch")
library(SynSigGen) # remotes::install_github(repo = "steverozen/SynSigGen", ref = "1.1.1-branch")

##################################################################
##                      Data preprocessing                      ##
##################################################################
# Get the real exposures from PCAWG assignments
real_exposures_sbs96 <- PCAWG7::exposure$PCAWG$SBS96
pcawg_sbs96_catalog <- PCAWG7::spectra$PCAWG$SBS96

# Only select samples that belong to the selected cancer types
cancer_types <- c(
  "Biliary-AdenoCA", "Breast-AdenoCA", "CNS-Medullo", 
  "ColoRect-AdenoCA", "Eso-AdenoCA", "Head-SCC",
  "Kidney-RCC", "Liver-HCC", "Lung-AdenoCA",
  "Lymph-BNHL", "Lymph-CLL", "Ovary-AdenoCA", 
  "Panc-AdenoCA", "Panc-Endocrine", "Prost-AdenoCA",
  "Skin-Melanoma", "Stomach-AdenoCA", "Uterus-AdenoCA"
)
indices_selected_types <- unlist(sapply(cancer_types, FUN = function(x) {
  grep(x, colnames(real_exposures_sbs96))
}))
real_exposures_sbs96 <- real_exposures_sbs96[, indices_selected_types]

# Exclude samples which have mutations less than 100
samples_less_than_100 <- names(which(colSums(pcawg_sbs96_catalog) < 100))
indices_less_than_100 <-
  which(colnames(real_exposures_sbs96) %in% samples_less_than_100)
real_exposures_sbs96 <-
  real_exposures_sbs96[, -indices_less_than_100, drop = FALSE]

# Exclude tumors which have possible artifact signatures
artifact_sigs <- cosmicsig::possible_artifacts()
indices_artifact_sigs <-
  which(rownames(real_exposures_sbs96) %in% artifact_sigs)
artifact_sigs_selected_types <-
  rownames(real_exposures_sbs96)[indices_artifact_sigs]
tumors_to_remove <- sapply(artifact_sigs_selected_types, FUN = function(x) {
  exposure <- real_exposures_sbs96[x, ]
  return(names(exposure[exposure > 0]))
})
tumors_to_remove <- unique(unlist(tumors_to_remove))
real_exposures_sbs96 <-
  real_exposures_sbs96[, !colnames(real_exposures_sbs96) %in% tumors_to_remove]

# Exclude signatures which have active tumors less than 10
sig_active_tumors <- apply(real_exposures_sbs96, MARGIN = 1, FUN = function(x) {
  return(length(x[x > 0]))
})
sigs_low_activity <-
  names(sort(sig_active_tumors[sig_active_tumors < 10], decreasing = TRUE))
samples_with_low_activity_sigs <- sapply(sigs_low_activity, FUN = function(x) {
  exposure <- real_exposures_sbs96[x, ]
  exposure_non_zero <- exposure[exposure > 0]
  return(names(exposure_non_zero))
})
samples_with_low_activity_sigs <- unique(unlist(samples_with_low_activity_sigs))
indices_low_activity_sigs <-
  which(colnames(real_exposures_sbs96) %in% samples_with_low_activity_sigs)
real_exposures_sbs96 <-
  real_exposures_sbs96[, -indices_low_activity_sigs, drop = FALSE]

real_exposures_sbs96 <- remove_zero_activity_sigs(real_exposures_sbs96)

##################################################################
##   Calculate number of synthetic tumors in each cancer type   ##
##################################################################

msi_sample_indices_selected_types <-
  unlist(sapply(pcawg_msi_tumor_ids, FUN = function(x) {
    grep(x, colnames(real_exposures_sbs96))
  }))
msi_sample_ids <- names(msi_sample_indices_selected_types)
length(msi_sample_ids) # 23

pole_sample_indices_selected_types <-
  unlist(sapply(pcawg_pole_tumor_ids, FUN = function(x) {
    grep(x, colnames(real_exposures_sbs96))
  }))
pole_sample_ids <- names(pole_sample_indices_selected_types)

# There are no POLE samples in real_exposures_sbs96
length(pole_sample_ids) # 0

real_exposures_sbs96_no_msi <-
  real_exposures_sbs96[, -msi_sample_indices_selected_types, drop = FALSE]
real_exposures_sbs96_no_msi <-
  remove_zero_activity_sigs(real_exposures_sbs96_no_msi)
real_exposures_sbs96_msi <-
  real_exposures_sbs96[, msi_sample_indices_selected_types, drop = FALSE]
real_exposures_sbs96_msi <- remove_zero_activity_sigs(real_exposures_sbs96_msi)

num_samples_total <- calculate_num_samples(real_exposures_sbs96)
num_samples_msi <- calculate_num_samples(real_exposures_sbs96_msi)
cancer_types_msi <- names(num_samples_msi)
num_samples_no_msi <- calculate_num_samples(real_exposures_sbs96_no_msi)

# Only generate 60 synthetic tumors for the selected cancer types 
# (total 1080). Scale the original number of tumors in each cancer type
# in real exposure accordingly
scale_factors <- 60 / num_samples_total

# Calculate the number of MSI-H synthetic tumors in each cancer type
num_samples_msi_scaled <-
  sapply(cancer_types_msi, FUN = function(x) {
    scaled_number <- ceiling(scale_factors[x] * num_samples_msi[x])
    names(scaled_number) <- NULL
    return(scaled_number)
  })

# Calculate the number of non MSI-H synthetic tumors in each cancer type
num_samples_no_msi_scaled <- rep(60, length(cancer_types))
names(num_samples_no_msi_scaled) <- cancer_types
for (i in cancer_types_msi) {
  num_samples_no_msi_scaled[i] <-
    num_samples_no_msi_scaled[i] - num_samples_msi_scaled[i]
}

# Make sure the total number of synthetic tumors is 1080
sum(num_samples_msi_scaled) + sum(num_samples_no_msi_scaled)

##################################################################
##                 Generation of synthetic data                 ##
##################################################################

output_dir_sbs96_no_msi <- "./SBS_set2/input/raw/PCAWG.SBS96.syn.exposures.no.msi"
output_dir_sbs96_msi <- "./SBS_set2/input/raw/PCAWG.SBS96.syn.exposures.msi"
output_dir_sbs96 <- "./SBS_set2/input/raw/PCAWG.SBS96.syn.exposures.no.noise"
output_dir_sbs96_nb_size_30 <-
  "./SBS_set2/input/raw/PCAWG.SBS96.syn.exposures.noisy.neg.binom.size.30"

distribution <- "neg.binom"
sample_prefix_name <- "SP.Syn."
mutation_type <- "SBS96"
seed <- 658220
input_sigs_sbs96 <- cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS96

sig_params_sbs96_selected_types <-
  SynSigGen:::GetSynSigParamsFromExposures(
    exposures = real_exposures_sbs96,
    distribution = distribution,
    sig.params = SynSigGen::signature.params$SBS96
  )

# Generate non MSI-H synthetic tumors
synthetic_tumors_sbs96_no_msi <-
  SynSigGen::GenerateSyntheticTumors(
    seed = seed,
    dir = output_dir_sbs96_no_msi,
    cancer.types = cancer_types,
    samples.per.cancer.type = num_samples_no_msi_scaled,
    input.sigs = input_sigs_sbs96,
    real.exposures = real_exposures_sbs96_no_msi,
    distribution = distribution,
    sample.prefix.name = sample_prefix_name,
    sig.params = sig_params_sbs96_selected_types
  )
unlink(output_dir_sbs96_no_msi, recursive = TRUE)
syn_exposures_sbs96_no_msi <-
  synthetic_tumors_sbs96_no_msi$ground.truth.exposures

# Generate MSI-H synthetic tumors
synthetic_tumors_sbs96_msi <-
  generate_subtype_syn_tumors(
    seed = seed,
    dir = output_dir_sbs96_msi,
    cancer_types = cancer_types_msi,
    samples_per_caner_type = num_samples_msi_scaled,
    input_sigs = input_sigs_sbs96,
    real_exposure = real_exposures_sbs96_msi,
    distribution = distribution,
    sample_prefix_name = sample_prefix_name,
    tumor_marker_name = "MSI-H",
    sig_params = sig_params_sbs96_selected_types
  )
unlink(output_dir_sbs96_msi, recursive = TRUE)
syn_exposures_sbs96_msi <-
  synthetic_tumors_sbs96_msi$ground.truth.exposures

# Combine the non MSI-H and MSI-H synthetic exposures in each cancer type
synthetic_exposures_sbs96 <-
  combine_exposure(
    syn_exposures_sbs96_no_msi,
    syn_exposures_sbs96_msi
  )

# Generate the combined synthetic tumors
write_sig_params(
  dir = output_dir_sbs96,
  real_exposure = real_exposures_sbs96,
  synthetic_exposure = synthetic_exposures_sbs96,
  cancer_types = cancer_types,
  distribution = distribution,
  sig_params = sig_params_sbs96_selected_types,
  sample_prefix_name = sample_prefix_name,
  mutation_type = mutation_type
)

catalog <- SynSigGen::CreateAndWriteCatalog(
  sigs = input_sigs_sbs96,
  exp = synthetic_exposures_sbs96,
  my.dir = output_dir_sbs96,
  overwrite = TRUE,
  extra.file.suffix = mutation_type
)

# Add noise to the synthetic tumors
sbs96_noisy_tumors_size_30 <-
  SynSigGen::GenerateNoisyTumors(
    seed = seed,
    dir = output_dir_sbs96_nb_size_30,
    input.exposure = synthetic_exposures_sbs96,
    signatures = input_sigs_sbs96,
    n.binom.size = 30
  )
noisy_exposures_size_30_sbs96 <- sbs96_noisy_tumors_size_30$exposures

#################################################################
##                   Plot data distributions                   ##
#################################################################

data_distribution_file <-
  "./SBS_set2/input/SBS_syn_data_distribution.pdf"
grDevices::pdf(
  file = data_distribution_file,
  width = 8.2677, height = 11.6929, onefile = TRUE
)
par(mfrow = c(4, 3), oma = c(0, 0, 1, 0))
plot_exposure_distribution(
  real_exposure = real_exposures_sbs96,
  synthetic_exposure = synthetic_exposures_sbs96,
  noisy_exposure = noisy_exposures_size_30_sbs96,
  size = 30,
  distribution = distribution,
  sig_params = sig_params_sbs96_selected_types,
  sample_prefix_name = sample_prefix_name
)
grDevices::dev.off()
