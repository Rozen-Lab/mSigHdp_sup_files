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
real_exposures_id <- PCAWG7::exposure$PCAWG$ID
pcawg_id_catalog <- PCAWG7::spectra$PCAWG$ID

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
  grep(x, colnames(real_exposures_id))
}))
real_exposures_id <- real_exposures_id[, indices_selected_types]

# Exclude samples which have mutations less than 100
samples_less_than_100 <- names(which(colSums(pcawg_id_catalog) < 100))
indices_less_than_100 <-
  which(colnames(real_exposures_id) %in% samples_less_than_100)
real_exposures_id <-
  real_exposures_id[, -indices_less_than_100, drop = FALSE]

# Exclude signatures which have active tumors less than 10
sig_active_tumors <- apply(real_exposures_id, MARGIN = 1, FUN = function(x) {
  return(length(x[x > 0]))
})
sigs_low_activity <-
  names(sort(sig_active_tumors[sig_active_tumors < 10], decreasing = TRUE))
samples_with_low_activity_sigs <- sapply(sigs_low_activity, FUN = function(x) {
  exposure <- real_exposures_id[x, ]
  exposure_non_zero <- exposure[exposure > 0]
  return(names(exposure_non_zero))
})
samples_with_low_activity_sigs <- unique(unlist(samples_with_low_activity_sigs))
indices_low_activity_sigs <-
  which(colnames(real_exposures_id) %in% samples_with_low_activity_sigs)
real_exposures_id <-
  real_exposures_id[, -indices_low_activity_sigs, drop = FALSE]

real_exposures_id <- remove_zero_activity_sigs(real_exposures_id)

##################################################################
##   Calculate number of synthetic tumors in each cancer type   ##
##################################################################

msi_sample_indices_selected_types <-
  unlist(sapply(pcawg_msi_tumor_ids, FUN = function(x) {
    grep(x, colnames(real_exposures_id))
  }))
msi_sample_ids <- names(msi_sample_indices_selected_types)
length(msi_sample_ids) # 37

pole_sample_indices_selected_types <-
  unlist(sapply(pcawg_pole_tumor_ids, FUN = function(x) {
    grep(x, colnames(real_exposures_id))
  }))
pole_sample_ids <- names(pole_sample_indices_selected_types)
length(pole_sample_ids) # 9

real_exposures_id_no_msi_pole <-
  real_exposures_id[, -c(
    msi_sample_indices_selected_types,
    pole_sample_indices_selected_types
  ), drop = FALSE]
real_exposures_id_no_msi_pole <-
  remove_zero_activity_sigs(real_exposures_id_no_msi_pole)
real_exposures_id_msi <-
  real_exposures_id[, msi_sample_indices_selected_types, drop = FALSE]
real_exposures_id_msi <- remove_zero_activity_sigs(real_exposures_id_msi)
real_exposures_id_pole <-
  real_exposures_id[, pole_sample_indices_selected_types, drop = FALSE]
real_exposures_id_pole <- remove_zero_activity_sigs(real_exposures_id_pole)

num_samples_total <- calculate_num_samples(real_exposures_id)
sum(num_samples_total) # 2019

num_samples_msi <- calculate_num_samples(real_exposures_id_msi)
cancer_types_msi <- names(num_samples_msi)

num_samples_pole <- calculate_num_samples(real_exposures_id_pole)
cancer_types_pole <- names(num_samples_pole)

num_samples_no_msi_pole <- calculate_num_samples(real_exposures_id_no_msi_pole)

# The number of synthetic tumors will be the same as the one in real exposure
sum(num_samples_msi) + sum(num_samples_pole) +
  sum(num_samples_no_msi_pole) # 2019

##################################################################
##                 Generation of synthetic data                 ##
##################################################################

output_dir_id_no_msi_pole <- "./indel_set2/input/raw/PCAWG.ID.syn.exposures.no.msi.pole"
output_dir_id_msi <- "./indel_set2/input/raw/PCAWG.ID.syn.exposures.msi"
output_dir_id_pole <- "./indel_set2/input/raw/PCAWG.ID.syn.exposures.pole"
output_dir_id <- "./indel_set2/input/raw/PCAWG.ID.syn.exposures.no.noise"
output_dir_id_nb_size_10 <-
  "./indel_set2/input/raw/PCAWG.ID.syn.exposures.noisy.neg.binom.size.10"

distribution <- "neg.binom"
sample_prefix_name <- "SP.Syn."
mutation_type <- "ID"
seed <- 658220
input_sigs_id <- cosmicsig::COSMIC_v3.2$signature$GRCh37$ID

sig_params_id_selected_types <-
  SynSigGen:::GetSynSigParamsFromExposures(
    exposures = real_exposures_id,
    distribution = distribution,
    sig.params = SynSigGen::signature.params$ID
  )

# Generate synthetic tumors that are not MSI-H or POLE-mutated
syn_tumors_id_no_msi_pole <-
  SynSigGen::GenerateSyntheticTumors(
    seed = seed,
    dir = output_dir_id_no_msi_pole,
    cancer.types = cancer_types,
    samples.per.cancer.type = num_samples_no_msi_pole,
    input.sigs = input_sigs_id,
    real.exposures = real_exposures_id_no_msi_pole,
    distribution = distribution,
    sample.prefix.name = sample_prefix_name,
    sig.params = sig_params_id_selected_types
  )
unlink(output_dir_id_no_msi_pole, recursive = TRUE)
syn_exposures_id_no_msi_pole <-
  syn_tumors_id_no_msi_pole$ground.truth.exposures

# Generate MSI-H synthetic tumors
synthetic_tumors_id_msi <-
  generate_subtype_syn_tumors(
    seed = seed,
    dir = output_dir_id_msi,
    cancer_types = cancer_types_msi,
    samples_per_caner_type = num_samples_msi,
    input_sigs = input_sigs_id,
    real_exposure = real_exposures_id_msi,
    distribution = distribution,
    sample_prefix_name = sample_prefix_name,
    tumor_marker_name = "MSI-H",
    sig_params = sig_params_id_selected_types
  )
unlink(output_dir_id_msi, recursive = TRUE)
syn_exposures_id_msi <-
  synthetic_tumors_id_msi$ground.truth.exposures

# Generate POLE-mutated synthetic tumors
synthetic_tumors_id_pole <-
  generate_subtype_syn_tumors(
    seed = seed,
    dir = output_dir_id_pole,
    cancer_types = cancer_types_pole,
    samples_per_caner_type = num_samples_pole,
    input_sigs = input_sigs_id,
    real_exposure = real_exposures_id_pole,
    distribution = distribution,
    sample_prefix_name = sample_prefix_name,
    tumor_marker_name = "POLE",
    sig_params = sig_params_id_selected_types
  )
unlink(output_dir_id_pole, recursive = TRUE)
syn_exposures_id_pole <-
  synthetic_tumors_id_pole$ground.truth.exposures

# Combine the non MSI-H non POLE, MSI-H and POLE synthetic exposures in each
# cancer type
synthetic_exposures_id <-
  combine_exposure(
    syn_exposures_id_no_msi_pole,
    syn_exposures_id_pole,
    syn_exposures_id_msi
  )

# Generate the combined synthetic tumors
write_sig_params(
  dir = output_dir_id,
  real_exposure = real_exposures_id,
  synthetic_exposure = synthetic_exposures_id,
  cancer_types = cancer_types,
  distribution = distribution,
  sig_params = sig_params_id_selected_types,
  sample_prefix_name = sample_prefix_name,
  mutation_type = mutation_type
)

catalog <- SynSigGen::CreateAndWriteCatalog(
  sigs = input_sigs_id,
  exp = synthetic_exposures_id,
  my.dir = output_dir_id,
  overwrite = TRUE,
  extra.file.suffix = mutation_type
)

# Add noise to the synthetic tumors
id_noisy_tumors_size_10 <-
  SynSigGen::GenerateNoisyTumors(
    seed = seed,
    dir = output_dir_id_nb_size_10,
    input.exposure = synthetic_exposures_id,
    signatures = input_sigs_id,
    n.binom.size = 10
  )
noisy_exposures_size_10_id <- id_noisy_tumors_size_10$exposures

#################################################################
##                   Plot data distributions                   ##
#################################################################

data_distribution_file <-
  "./indel_set2/input/indel_syn_data_distributions.pdf"
grDevices::pdf(
  file = data_distribution_file,
  width = 8.2677, height = 11.6929, onefile = TRUE
)
par(mfrow = c(4, 3), oma = c(0, 0, 1, 0))
plot_exposure_distribution(
  real_exposure = real_exposures_id,
  synthetic_exposure = synthetic_exposures_id,
  noisy_exposure = noisy_exposures_size_10_id,
  size = 10,
  distribution = distribution,
  sig_params = sig_params_id_selected_types,
  sample_prefix_name = sample_prefix_name
)
grDevices::dev.off()
