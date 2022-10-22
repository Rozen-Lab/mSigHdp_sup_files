# This file should be sourced from
# SBS_set2/code/1_data_generation.R
# indel_set2/code/1_data_generation.R
###########################################################################

# Source utility functions shared by "data_gen_utils_set1.R" 
# and "data_gen_utils_set2.R"
basedir <- "mSigHdp_sup_files_test" 
if (basename(getwd()) != basedir) {
  stop("Please run from top level directory, ", basedir)
}
source("./common_code/data_gen_utils.R")


#' Plot the distributions of several exposures
#'
#' @param real_exposure A matrix of real exposures.
#'
#' @param synthetic_exposure A matrix of synthetic exposures.
#'
#' @param noisy_exposure A matrix of synthetic exposures with noise
#'   added.
#'
#' @param size Negative binomial size parameter used for adding noise to
#'   \code{synthetic_exposure}.
#'
#' @param distribution Probability distribution used to model exposures due to
#'   one mutational signature. Can be \code{neg.binom} which stands for negative
#'   binomial distribution.
#'
#' @param sig_params Empirical signature parameters generated using real
#'   exposures irrespective of their cancer types. If there is only \strong{one}
#'   tumor having a signature in a cancer type in \code{real_exposure}, we
#'   cannot fit the \code{distribution} to only one data point. Instead, we will
#'   use the empirical parameter \code{size} from \code{sig_params}. The
#'   mutation of the only one tumor will be used as parameter \code{mu} in
#'   negative binomial distribution.
#'
#' @param sample_prefix_name Prefix name to add to the synthetic tumors.
plot_exposure_distribution <-
  function(real_exposure, synthetic_exposure,
           noisy_exposure, size,
           distribution,
           sig_params, sample_prefix_name) {
    real_exposure_info <-
      get_exposure_info(
        exposure = real_exposure,
        distribution = distribution,
        sig_params = sig_params
      )
    synthetic_exposure_info <-
      get_exposure_info(
        exposure = synthetic_exposure,
        distribution = distribution,
        sig_params = sig_params
      )
    noisy_exposure_info <-
      get_exposure_info(
        exposure = noisy_exposure,
        distribution = distribution,
        sig_params = sig_params
      )
    
    cancer_types <- names(real_exposure_info$exposures)
    
    retval <- sapply(cancer_types, FUN = function(x) {
      one_caner_type <- x
      one_caner_type_with_prefix <- paste0(sample_prefix_name, one_caner_type)
      real_exposure_one_type <- real_exposure_info$exposures[[one_caner_type]]
      synthetic_exposure_one_type <-
        synthetic_exposure_info$exposures[[one_caner_type_with_prefix]]
      noisy_exposure_one_type <-
        noisy_exposure_info$exposures[[one_caner_type_with_prefix]]
      
      params_real_one_type <- real_exposure_info$sig_params[[one_caner_type]]
      params_synthetic_one_type <-
        synthetic_exposure_info$sig_params[[one_caner_type_with_prefix]]
      params_noisy_one_type <-
        noisy_exposure_info$sig_params[[one_caner_type_with_prefix]]
      
      plot_exposure_dist_one_type(
        real_exposure = real_exposure_one_type,
        synthetic_exposure = synthetic_exposure_one_type,
        params_real = params_real_one_type,
        params_synthetic = params_synthetic_one_type,
        noisy_exposure = noisy_exposure_one_type,
        params_noisy = params_noisy_one_type,
        size = size,
        cancer_type = one_caner_type
      )
    })
  }


#' Plot the distributions of exposures attributed to mutational signatures in
#' one cancer type
#'
#' @param real_exposure A matrix of real exposures.
#'
#' @param synthetic_exposure A matrix of synthetic exposures.
#'
#' @param params_real Signature parameters of \code{real_exposure}.
#'
#' @param params_synthetic Signature parameters of \code{synthetic_exposure}.
#'
#' @param noisy_exposure A matrix of synthetic exposures with noise added.
#'
#' @param params_noisy Signature parameters of \code{noisy_exposure}.
#'
#' @param size Negative binomial size parameter used for adding noise to
#'   \code{synthetic_exposure}.
#'
#' @param cancer_type A character string denoting one cancer type. See
#'   \code{PCAWG7::CancerTypes()} for examples.
#'
plot_exposure_dist_one_type <-
  function(real_exposure, synthetic_exposure,
           params_real, params_synthetic,
           noisy_exposure, params_noisy, size,
           cancer_type) {
    real_exposure <- remove_zero_activity_sigs(real_exposure)
    synthetic_exposure <- remove_zero_activity_sigs(synthetic_exposure)
    noisy_exposure <- remove_zero_activity_sigs(noisy_exposure)
    
    for (i in rownames(synthetic_exposure)) {
      real_exposure_non_zero <- real_exposure[i, which(real_exposure[i, ] > 0)]
      synthetic_exposure_non_zero <-
        synthetic_exposure[i, which(synthetic_exposure[i, ] > 0)]
      if (i %in% rownames(noisy_exposure)) {
        noisy_exposure_non_zero <-
          noisy_exposure[i, which(noisy_exposure[i, ] > 0)]
      } else {
        noisy_exposure_non_zero <- NULL
      }
      
      draw_histogram(
        counts = real_exposure_non_zero,
        title = paste0(cancer_type, ".", i, ".real.exposure"),
        cex_main = 0.9,
        params = params_real,
        sig_id = i
      )
      
      draw_histogram(
        counts = synthetic_exposure_non_zero,
        title = paste0(cancer_type, ".", i, ".synthetic.exposure"),
        cex_main = 0.9,
        params = params_synthetic,
        sig_id = i
      )
      
      if (!is.null(noisy_exposure_non_zero)) {
        draw_histogram(
          counts = noisy_exposure_non_zero,
          title = paste0(cancer_type, ".", i, ".noisy.exposure"),
          cex_main = 0.85,
          params = params_noisy,
          sig_id = i,
          size = size
        )
      } else {
        plot.new()
      }
      
      if (length(real_exposure_non_zero) > 1) {
        draw_density(
          counts = real_exposure_non_zero,
          title = paste0(cancer_type, ".", i, ".real.exposure"),
          cex_main = 0.9, params = params_real, sig_id = i,
        )
      } else {
        plot.new()
      }
      
      # Plot several densities in one graph
      if (length(real_exposure_non_zero) > 1) {
        draw_two_density(
          counts1 = real_exposure_non_zero,
          counts2 = synthetic_exposure_non_zero,
          title = paste0(
            cancer_type, ".", i, ".synthetic.exposure"
          ),
          cex_main = 0.9,
          params = params_synthetic,
          sig_id = i,
          legend = c("real.exposure", "synthetic.exposure")
        )
      } else {
        plot.new()
      }
      
      if (length(real_exposure_non_zero) > 1) {
        draw_two_density(
          counts1 = real_exposure_non_zero,
          counts2 = noisy_exposure_non_zero,
          title = paste0(
            cancer_type, ".", i, ".noisy.exposure"
          ),
          cex_main = 0.85,
          params = params_noisy,
          sig_id = i,
          legend = c("real.exposure", "noisy.exposure"),
          size = size
        )
      } else {
        plot.new()
      }
    }
  }

#' Create multiple box plots showing the scaled Euclidean distances between
#' different groups
#'
#' @param distance_df A data frame which has the scaled Euclidean distances
#'   information from multiple groups.
#'
#' @param data_type Mutation type of the data which generated \code{distance_df}
#'   (e.g. "SBS", "indel").
#'
#' @param ylim Y-axis limit of the box plot.
#'
#' @return A list of multiple ggplot objects.
#'
create_boxplots <- function(distance_df, data_type, ylim) {
  cancer_types <- unique(distance_df$cancer.type)
  # First create box plot for all cancer types
  plot_objects <-
    one_boxplot(
      distance_df = distance_df,
      title = paste0("All 18 cancer types (", data_type, ")"),
      ylim = ylim
    )
  
  # Create box plots for individual cancer type
  for (cancer_type in cancer_types) {
    distance_df_one_type <-
      distance_df[distance_df$cancer.type == cancer_type, ]
    tmp <-
      one_boxplot(
        distance_df = distance_df_one_type,
        title = paste0(cancer_type, " (", data_type, ")"),
        ylim = ylim
      )
    plot_objects <- c(plot_objects, tmp)
  }
  
  return(plot_objects)
}
