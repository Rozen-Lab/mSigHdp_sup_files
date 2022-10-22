# This file should be sourced from
# common_code/data_gen_utils_set1.R
# common_code/data_gen_utils_set2.R
# common_code/noise_selection/code/SBS_noise_selection.R
# common_code/noise_selection/code/indel_noise_selection.R
###########################################################################

#' Remove signatures which have zero activity from an exposure matrix
#'
#' @param exposure Exposures as a numerical matrix (or data.frame) with
#'   signatures in rows and samples in columns. Rownames are taken as the
#'   signature names and column names are taken as the sample IDs.
#'
#' @return A matrix (or data.frame) with zero activity signatures removed.
#'
remove_zero_activity_sigs <- function(exposure) {
  return(exposure[rowSums(exposure) > 0, , drop = FALSE])
}

#' Calculate number of samples for each cancer type in an exposure matrix
#'
#' @param exposure Exposures as a numerical matrix (or data.frame) with
#'   signatures in rows and samples in columns. Rownames are taken as the
#'   signature names and column names are taken as the sample IDs.
#'
#' @importFrom PCAWG7 SplitPCAWGMatrixByTumorType
#'
#' @return A vector of integers showing the number of samples for each cancer
#'   type from \code{exposure}.
#'
calculate_num_samples <- function(exposure) {
  exposures_cancer_types <- PCAWG7::SplitPCAWGMatrixByTumorType(exposure)
  return(sapply(exposures_cancer_types, FUN = ncol))
}

#' Generate synthetic tumors from a subtype of exposure matrix
#'
#' @param seed A random seed to use.
#'
#' @param dir The directory in which to put the output; will be created if
#'   necessary.
#'
#' @param cancer_types A vector of character strings denoting different cancer
#'   types. This function will search \code{real_exposure} for exposures from
#'   tumors matching these strings. See \code{PCAWG7::CancerTypes()} for
#'   example.
#'
#' @param samples_per_caner_type Number of synthetic tumors to create for each
#'   cancer type. If it is \strong{one} number, then generate the \strong{same}
#'   number of synthetic tumors for each \code{cancer_types}. Or if it is a
#'   \strong{vector} of numbers, then generate synthetic tumors for each
#'   \code{cancer_types} accordingly to the number specified in the vector. The
#'   length and order of \code{samples_per_caner_type} should match that in
#'   \code{cancer_types}.
#'
#' @param input_sigs A matrix of signatures.
#'
#' @param real_exposure A matrix of real exposures.
#'
#' @param distribution Probability distribution used to model exposures due to
#'   one mutational signature. Can be \code{neg.binom} which stands for negative
#'   binomial distribution.
#'
#' @param sample_prefix_name Prefix name to add to the synthetic tumors.
#'
#' @param tumor_marker_name Tumor marker name to add to the synthetic tumors.
#'   E.g. "MSI-H", "POLE".
#'
#' @param sig_params Empirical signature parameters generated using real
#'   exposures irrespective of their cancer types. If there is only \strong{one}
#'   tumor having a signature in a cancer type in \code{real_exposure}, we
#'   cannot fit the \code{distribution} to only one data point. Instead, we will
#'   use the empirical parameter \code{size} from \code{sig_params}. The
#'   mutation of the only one tumor will be used as parameter \code{mu} in
#'   negative binomial distribution.
#'
#' @importFrom SynSigGen GetSynSigParamsFromExposures GenerateListOfSigParams
#'   GenerateSyntheticTumorsFromSigParams
#'
#' @return A list of three elements that comprise the synthetic data:
#' \enumerate{
#'  \item \code{ground.truth.catalog}: Spectra catalog with rows denoting
#'  mutation types and columns denoting sample names.
#'
#'  \item \code{ground.truth.signatures}: Signatures active in
#'  \code{ground.truth.catalog}.
#'
#'  \item \code{ground.truth.exposures}: Exposures of
#'  \code{ground.truth.signatures} in \code{ground.truth.catalog}.
#' }
#'
generate_subtype_syn_tumors <-
  function(seed, dir, cancer_types, samples_per_caner_type, input_sigs,
           real_exposure, distribution, sample_prefix_name, tumor_marker_name,
           sig_params) {
    subtype_sig_params_all <-
      SynSigGen::GetSynSigParamsFromExposures(
        exposures = real_exposure,
        distribution = distribution,
        cancer.type = tumor_marker_name,
        sig.params = sig_params
      )

    subtype_sig_params <-
      SynSigGen::GenerateListOfSigParams(
        real.exposures = real_exposure,
        cancer.types = cancer_types,
        distribution = distribution,
        sig.params = subtype_sig_params_all
      )

    synthetic_tumors <-
      SynSigGen::GenerateSyntheticTumorsFromSigParams(
        seed = seed,
        dir = dir,
        cancer.types = cancer_types,
        samples.per.cancer.type = samples_per_caner_type,
        input.sigs = input_sigs,
        sig.params = subtype_sig_params,
        distribution = distribution,
        sample.prefix.name = sample_prefix_name,
        tumor.marker.name = tumor_marker_name,
      )
    return(synthetic_tumors)
  }

#' Combine several exposure matrix according to cancer types
#'
#' @param exposure1 Exposures as a numerical matrix (or data.frame) with
#'   signatures in rows and samples in columns. Rownames are taken as the
#'   signature names and column names are taken as the sample IDs.
#'
#' @param exposure2 Exposures as a numerical matrix (or data.frame) with
#'   signatures in rows and samples in columns. Rownames are taken as the
#'   signature names and column names are taken as the sample IDs.
#'
#' @param exposure3 Exposures as a numerical matrix (or data.frame) with
#'   signatures in rows and samples in columns. Rownames are taken as the
#'   signature names and column names are taken as the sample IDs.
#'
#' @section Note
#' \code{exposure1} should have more cancer types compared to \code{exposure2}
#' or \code{exposure3}.
#'
#' @importFrom PCAWG7 SplitPCAWGMatrixByTumorType
#'
#' @importFrom SynSigGen MergeExposures
#'
#' @return A combined exposure matrix according to cancer type.
#'
combine_exposure <- function(exposure1, exposure2, exposure3 = NULL) {
  exposure1_cancer_types <- PCAWG7::SplitPCAWGMatrixByTumorType(exposure1)
  exposure2_cancer_types <- PCAWG7::SplitPCAWGMatrixByTumorType(exposure2)
  if (!is.null(exposure3)) {
    exposure3_cancer_types <- PCAWG7::SplitPCAWGMatrixByTumorType(exposure3)
  }

  exposure_list <-
    lapply(names(exposure1_cancer_types), FUN = function(x) {
      cancer_type <- x
      exposure1_cancer_type <- exposure1_cancer_types[[cancer_type]]
      exposure2_cancer_type <- exposure2_cancer_types[[cancer_type]]

      if (!is.null(exposure3)) {
        exposure3_cancer_type <- exposure3_cancer_types[[cancer_type]]
        list_of_exposure <- list(
          exposure1_cancer_type, exposure2_cancer_type,
          exposure3_cancer_type
        )
      } else {
        list_of_exposure <- list(exposure1_cancer_type, exposure2_cancer_type)
      }

      list_of_exposure2 <- Filter(Negate(is.null), list_of_exposure)
      merged_exposure <-
        SynSigGen::MergeExposures(list.of.exposures = list_of_exposure2)
      return(merged_exposure)
    })
  new_exposure <- SynSigGen::MergeExposures(list.of.exposures = exposure_list)
  return(new_exposure)
}

#' Write the signature parameters of real and synthetic exposures
#'
#' @inheritParams generate_subtype_syn_tumors
#'
#' @param synthetic_exposure A matrix of synthetic exposures.
#'
#' @param mutation_type Mutation type of the data, e.g. "SBS96", "ID".
#'
#' @importFrom PCAWG7 SplitPCAWGMatrixByTumorType
#'
#' @importFrom SynSigGen GenerateListOfSigParams WriteSynSigParams
#'   GetSynSigParamsFromExposures
#'
write_sig_params <-
  function(dir, real_exposure, synthetic_exposure, cancer_types, distribution,
           sig_params, sample_prefix_name, mutation_type) {

    # Write parameters into files
    froot <- file.path(dir, "parameters")
    dir.create(path = froot, showWarnings = FALSE, recursive = TRUE)

    # Getting empirical estimates of key parameters describing exposures due to
    # signatures
    params <- SynSigGen::GenerateListOfSigParams(
      real.exposures = real_exposure,
      cancer.types = cancer_types,
      distribution = distribution,
      sig.params = sig_params
    )

    syn_exposure_cancer_types <-
      PCAWG7::SplitPCAWGMatrixByTumorType(synthetic_exposure)

    lapply(cancer_types, FUN = function(x) {
      one_cancer_type <- x
      parms <- params[[one_cancer_type]]

      parm_file <- file.path(froot, paste0(
        sample_prefix_name, one_cancer_type,
        ".params.", mutation_type, ".csv"
      ))
      cat("# Original paramaters\n", file = parm_file)
      # Suppress warning on column names on append
      suppressWarnings(
        SynSigGen::WriteSynSigParams(parms, parm_file,
          append = TRUE,
          col.names = NA
        )
      )
      index <- grep(
        pattern = one_cancer_type,
        x = names(syn_exposure_cancer_types)
      )
      syn_exp <- syn_exposure_cancer_types[[index]]

      # Sanity check; we regenerate the parameters from the synthetic exposures.
      check_params <-
        SynSigGen::GetSynSigParamsFromExposures(
          exposures = syn_exp,
          distribution = distribution,
          cancer.type = one_cancer_type,
          sig.params = sig_params
        )

      # check_params should be similar to parms
      cat("# Parameters derived from synthetic exposures\n",
        file = parm_file, append = TRUE
      )

      missing_sig_names <- setdiff(colnames(parms), colnames(check_params))
      if (length(missing_sig_names) > 0) {
        check_param2 <- matrix(NA, nrow = dim(parms)[1], ncol = dim(parms)[2])
        dimnames(check_param2) <- dimnames(parms)
        check_param2[, colnames(check_params)] <- check_params
        check_params <- check_param2
      } else {
        check_params <- check_params[, colnames(parms), drop = FALSE]
      }

      suppressWarnings(
        SynSigGen::WriteSynSigParams(check_params, parm_file,
          append = TRUE, col.names = NA
        )
      )

      if (length(missing_sig_names) > 0) {
        cat("# Some signatures not represented in the synthetic data:\n",
          file = parm_file, append = TRUE
        )
        cat("#", missing_sig_names, "\n", file = parm_file, append = TRUE)
      }

      cat("# Difference between original parameters and parameters",
        "derived from synthetic exposures\n",
        file = parm_file, append = TRUE
      )
      SynSigGen::WriteSynSigParams(parms - check_params, parm_file,
        append = TRUE, col.names = NA
      )
    })
    invisible(TRUE)
  }

#' Get the individual cancer type exposure and parameters from a large exposure
#' matrix
#'
#' @inheritParams generate_subtype_syn_tumors
#'
#' @importFrom PCAWG7 SplitPCAWGMatrixByTumorType
#'
#' @importFrom SynSigGen GenerateListOfSigParams
#'
#' @return A list of two elements:
#' \enumerate{
#'  \item \code{exposures}: A list of exposures for tumors in each cancer type.
#'  \item \code{sig_params}: A list of signature parameters for tumors in each
#'  cancer type.
#' }
get_exposure_info <- function(exposure, distribution, sig_params) {
  exposure_cancer_types <- PCAWG7::SplitPCAWGMatrixByTumorType(exposure)
  cancer_types <- names(exposure_cancer_types)
  sig_params_cancer_types <-
    SynSigGen::GenerateListOfSigParams(
      real.exposures = exposure,
      cancer.types = cancer_types,
      distribution = distribution,
      sig.params = sig_params
    )
  return(list(
    exposures = exposure_cancer_types,
    sig_params = sig_params_cancer_types
  ))
}

# Sample IDs of POLE-mutated tumors in the paper by Alexandrov, Kim,
# Haradhvala, Huang et al., 'The repertoire of Mutational Signatures in Human
# Cancer'. https://doi.org/10.1038/s41586-020-1943-3
pcawg_pole_tumor_ids <-
  c(
    "ColoRect-AdenoCA::SP81312", "ColoRect-AdenoCA::SP22031",
    "ColoRect-AdenoCA::SP16886", "ColoRect-AdenoCA::SP19295",
    "ColoRect-AdenoCA::SP17905", "ColoRect-AdenoCA::SP21400",
    "ColoRect-AdenoCA::SP18946", "ColoRect-AdenoCA::SP80615",
    "Uterus-AdenoCA::SP92659"
  )

# Sample IDs of MSI-H tumors in the paper by Alexandrov, Kim, Haradhvala, Huang
# et al., 'The repertoire of Mutational Signatures in Human Cancer'.
# https://doi.org/10.1038/s41586-020-1943-3
pcawg_msi_tumor_ids <-
  c(
    "Biliary-AdenoCA::SP99325", "ColoRect-AdenoCA::SP18310",
    "ColoRect-AdenoCA::SP17172", "ColoRect-AdenoCA::SP96133",
    "ColoRect-AdenoCA::SP110242", "ColoRect-AdenoCA::SP21017",
    "ColoRect-AdenoCA::SP22383", "ColoRect-AdenoCA::SP19215",
    "ColoRect-AdenoCA::SP18121", "ColoRect-AdenoCA::SP96118",
    "Kidney-RCC::SP102897", "Liver-HCC::SP98845",
    "Liver-HCC::SP107012", "Lymph-BNHL::SP116697",
    "Ovary-AdenoCA::SP102133", "Panc-AdenoCA::SP125732",
    "Panc-AdenoCA::SP125746", "Panc-AdenoCA::SP125770",
    "Skin-Melanoma::SP124384", "Stomach-AdenoCA::SP135287",
    "Stomach-AdenoCA::SP85122", "Stomach-AdenoCA::SP84384",
    "Stomach-AdenoCA::SP84982", "Stomach-AdenoCA::SP84392",
    "Stomach-AdenoCA::SP84439", "Uterus-AdenoCA::SP95454",
    "Uterus-AdenoCA::SP94741", "Uterus-AdenoCA::SP90209",
    "Uterus-AdenoCA::SP90989", "Uterus-AdenoCA::SP93540",
    "Uterus-AdenoCA::SP93227", "Uterus-AdenoCA::SP94348",
    "Uterus-AdenoCA::SP94933", "Uterus-AdenoCA::SP92364",
    "Uterus-AdenoCA::SP94917", "Uterus-AdenoCA::SP92460",
    "Uterus-AdenoCA::SP89909", "Kidney-ChRCC::SP123975"
  )

#' Get the size parameter information
#'
#' @param size The size parameter used in negative binomial distribution for
#'   adding noise.
#'
#' @return A character string showing the size parameter information.
get_size_message <- function(size) {
  if (is.null(size)) {
    return("")
  } else {
    return(paste0("\n neg.binom.size = ", size))
  }
}

#' Plot the histogram of mutation counts attributed to one mutational signature
#'
#' @param counts A vector of mutation counts attributed to one signature across
#'   length(counts) samples.
#'
#' @param title The first line of title on top of the histogram.
#'
#' @param cex_main The magnification to be used for main titles relative to the
#'   current setting of cex.
#'
#' @param params A matrix of parameters of exposures attributed to mutational
#'   signatures. The rows are
#'   \describe{
#'      \item{prob}{The proportion of tumors with the signature.}
#'
#'      \item{size}{Dispersion parameter of negative binomial distribution
#'      fitted to the exposures.}
#'
#'      \item{mu}{Mean of negative binomial distribution fitted to the
#'      exposures.}
#'   }
#'   See \code{SynSigGen::signature.params$SBS96} for example.
#'
#' @param sig_id Signature ID, e.g. "SBS1".
#'
#' @param size The size parameter used in negative binomial distribution for
#'   adding noise.
#'
draw_histogram <-
  function(counts, title, cex_main, params, sig_id, size = NULL) {
    size_message <- get_size_message(size)
    
    hist(counts,
         breaks = 1000,
         main = paste0(
           title,
           "\n N = ", length(counts),
           " prob = ", round(params[1, sig_id], 4), " ",
           size_message,
           "\n mu = ", round(params[3, sig_id], 2),
           "\n size = ", round(params[2, sig_id], 2)
         ),
         cex.main = cex_main,
         xlab = "counts",
         probability = TRUE
    )
    
    if (length(counts) > 1) {
      lines(density(counts), col = "red")
    }
  }

#' Draw kernel density plot of mutation counts attributed to one mutational
#' signature
#'
#' @inheritParams draw_histogram
#'
#' @param num_sample Number of samples to be shown at the title.
#'
#' @param xlab Label for the x axis.
#'
draw_density <- function(counts, title, cex_main, params, sig_id,
                         num_sample = length(counts),
                         size = NULL, xlab = NULL) {
  size_message <- get_size_message(size)
  
  plot(density(counts),
       col = "red",
       main = paste0(
         title,
         "\n N = ", num_sample,
         " prob = ", round(params[1, sig_id], 4), " ",
         size_message,
         "\n mu = ", round(params[3, sig_id], 2),
         "\n size = ", round(params[2, sig_id], 2)
       ),
       cex.main = cex_main,
       xlab = xlab
  )
}

#' Draw two kernel density plots of mutation counts attributed to mutational
#' signature
#'
#' @inheritParams draw_density
#'
#' @param counts1 A vector of mutation counts attributed to one signature across
#'   \code{length(counts1)} samples.
#'
#' @param counts2 A vector of mutation counts attributed to one signature across
#'   \code{length(counts2)} samples.
#'
#' @param legend A character or expression vector of length = 1 to appear in the
#'   legend.
#'
draw_two_density <-
  function(counts1, counts2, title, cex_main, params, sig_id,
           legend, size = NULL,
           xlab = "Two-sample Kolmogorov-Smirnov test") {
    draw_density(
      counts = counts1, title = title, cex_main = cex_main,
      params = params, sig_id = sig_id, num_sample = length(counts2),
      size = size, xlab = xlab
    )
    
    if (length(counts2) > 1) {
      lines(density(counts2),
            col = "blue",
            xlim = c(0, max(counts1))
      )
    }
    
    retval <-
      suppressWarnings(ks.test(
        x = counts1,
        y = counts2
      ))
    p_value <- round(retval$p.value, 3)
    legend("topright",
           title = paste0("p-value = ", p_value),
           legend = legend,
           col = c("red", "blue"),
           fill = c("red", "blue"),
           border = "white", bty = "n"
    )
  }

#' Plot catalog to a PDF file
#'
#' @param catalog A catalog as defined in \code{\link{ICAMS}} with attributes
#'   added. Use \code{ICAMS::ReadCatalog} to read in a catalog.
#'
#' @param identifier Character string showing the identity of \code{catalog},
#'   e.g. "no noise", "more noise", "less noise".
#'
#' @param file The name of the PDF file to be produced.
#'
#' @param grid A logical value indicating whether to draw grid lines.
#'
#' @param upper A logical value indicating whether to draw horizontal lines and
#'   the names of major mutation class on top of graph.
#'
#' @param xlabels A logical value indicating whether to draw x axis labels.
#'
plot_catalog_to_pdf <- function(catalog, identifier, file,
                                grid = FALSE, upper = FALSE, xlabels = FALSE) {
  colnames(catalog) <- paste0(
    colnames(catalog), " (", identifier, ") ", "(counts = ",
    colSums(catalog), ")"
  )
  
  if (inherits(catalog, "SBS96Catalog")) {
    # Suppress plotting  x axis tick marks for SBS96Catalog
    old_par <- par(tck = 0)
    on.exit(par(old_par))
  }
  
  ICAMS::PlotCatalogToPdf(
    catalog = catalog,
    file = file,
    grid = grid,
    upper = upper,
    xlabels = xlabels
  )
}

#' Get cancer type information from sample IDs of tumors
#'
#' @param sample_ids Character vector of sample IDs from tumors (e.g.
#'   "Breast-AdenoCA::SP117975", "SP.Syn.Breast-AdenoCA::S.1").
#'
#' @param syn_tumor Logical value indicating whether \code{sample_ids} are from
#'   synthetic tumors.
#'
#' @return Character vector of cancer types of tumors.
#'
get_cancer_type <- function(sample_ids, syn_tumor = FALSE) {
  cancer_types <- sapply(sample_ids, FUN = function(sample_id) {
    tmp <- strsplit(sample_id, "::")[[1]][1]
    if (syn_tumor) {
      return(strsplit(tmp, ".", fixed = TRUE)[[1]][3])
    } else {
      return(tmp)
    }
  })
  return(cancer_types)
}

#' Calculate scaled Euclidean distances between spectra and reconstructed
#' spectra of tumors
#'
#' @param spectra Spectra catalog as a numerical matrix with rows denoting
#'   mutation types and columns denoting sample names.
#'
#' @param exposure Exposures as a numerical matrix (or data.frame) with
#'   signatures in rows and samples in columns. Rownames are taken as the
#'   signature names and column names are taken as the sample IDs.
#'
#' @param sigs A matrix of signatures.
#'
#' @param group Character string showing the group information for
#'   \code{spectra}.
#'
#' @param syn_tumor Logical value indicating whether \code{spectra} are from
#'   synthetic tumors.
#'
#' @return A data frame with row names being the sample names in \code{spectra}.
#'   The columns show Euclidean distances information between \code{spectra} and
#'   reconstructed spectra using \code{sigs} and \code{exposure}. The column
#'   names are:
#'   * euclidean: Euclidean distance
#'   * mutations: Mutations of \code{spectra}
#'   * scaled.euclidean: euclidean / mutations
#'   * group: Group information for \code{spectra}
#'   * cancer.type: Cancer type for each sample in \code{spectra}
#'
get_distance <- function(spectra, exposure, sigs, group, syn_tumor = FALSE) {
  distance <- mSigAct:::CalculateDistance(
    spectra = spectra,
    exposure = exposure,
    sigs = sigs
  )
  distance <- distance[, c("euclidean", "mutations", "scaled.euclidean")]
  distance$group <- group
  distance$cancer.type <-
    get_cancer_type(sample_ids = rownames(distance), syn_tumor = syn_tumor)
  return(distance)
}

#' Read in synthetic data from a directory
#'
#' @param dir Directory path to read in the synthetic data.
#'
#' @return A list of three elements that comprise the synthetic data:
#' \enumerate{
#'  \item \code{spectra}: Spectra catalog with rows denoting
#'  mutation types and columns denoting sample names.
#'
#'  \item \code{exposure}: Exposures of \code{sigs} in \code{spectra}.
#'
#'  \item \code{sigs}: Signatures active in \code{spectra}.
#' }
#'
get_syn_data_info <- function(dir) {
  syn_spectra <-
    ICAMS::ReadCatalog(file = file.path(
      dir,
      "ground.truth.syn.catalog.csv"
    ))
  syn_exposure <-
    mSigAct::ReadExposure(file = file.path(
      dir,
      "ground.truth.syn.exposures.csv"
    ))
  sigs <-
    ICAMS::ReadCatalog(file = file.path(
      dir,
      "ground.truth.syn.sigs.csv"
    ), catalog.type = "counts.signature")
  return(list(
    spectra = syn_spectra, exposure = syn_exposure,
    sigs = sigs
  ))
}

#' Add noise to noiseless synthetic data with different negative-binomial size
#' parameter
#'
#' @param seed A random seed to use.
#'
#' @param exposure Noiseless synthetic exposure.
#'
#' @param sigs A matrix of signatures.
#'
#' @param nb_sizes Character vector of numerical values specifying different
#'   negative-binomial size parameter.
#'
#' @return A list with the same length as \code{nb_sizes}. Each element in the
#'   list is a list which contains three elements that comprise the synthetic
#'   data with noise added:
#'   \enumerate{
#'    \item \code{spectra}: Spectra catalog with noise added.
#'
#'    \item \code{exposure}: Exposures of \code{sigs} in \code{spectra}.
#'
#'    \item \code{sigs}: Signatures active in \code{spectra}.
#' }
#'
generate_noisy_data <- function(seed, exposure, sigs, nb_sizes) {
  retval <- lapply(nb_sizes, FUN = function(size) {
    output_dir <- file.path(tempdir(), size)
    tmp <-
      SynSigGen::GenerateNoisyTumors(
        seed = seed,
        dir = output_dir,
        input.exposure = exposure,
        signatures = sigs,
        n.binom.size = size
      )
    # Delete temp dir
    unlink(x = output_dir, recursive = TRUE)
    noise_data <- list(
      spectra = tmp$spectra, exposure = tmp$exposures,
      sigs = sigs
    )
    return(noise_data)
  })
  names(retval) <- paste0("size_", nb_sizes)
  return(retval)
}

#' Calculate scaled Euclidean distances for a list of synthetic data
#'
#' @param list_of_syn_data A list of synthetic data.
#'
#' @return A list of data frames showing the scaled Euclidean distance
#'   information. See the return value of function \code{get_distance}
#'   for more details.
#'
get_multiple_syn_distances <- function(list_of_syn_data) {
  retval <-
    lapply(seq_len(length(list_of_syn_data)), FUN = function(index) {
      one_syn_data <- list_of_syn_data[[index]]
      group <- names(list_of_syn_data)[index]
      out <- get_distance(
        spectra = one_syn_data$spectra,
        exposure = one_syn_data$exposure,
        sigs = one_syn_data$sigs,
        group = group,
        syn_tumor = TRUE
      )
      return(out)
    })
  names(retval) <- names(list_of_syn_data)
  return(retval)
}

#' Create one box plot showing the scaled Euclidean distances between different
#' groups
#'
#' @param distance_df A data frame which has the scaled Euclidean distances
#'   information from multiple groups.
#'
#' @param title Title of the box plot.
#'
#' @param ylim Y-axis limit of the box plot.
#'
#' @return A list of one ggplot object.
#'
one_boxplot <- function(distance_df, title, ylim) {
  # Customize x-axis tick labels
  x_lables <- unique(distance_df$group)
  x_lables <-
    gsub(pattern = "real", replacement = "Actual\ndata", x = x_lables)
  x_lables <-
    gsub(pattern = "size_", replacement = "Size=", x = x_lables)
  
  plot_object <-
    ggpubr::ggboxplot(
      data = distance_df,
      x = "group", y = "scaled.euclidean", color = "group"
    ) +
    ggpubr::stat_compare_means(
      label = "p.format", ref.group = "real", method = "wilcox.test",
      method.args = list(alternative = "two.sided")
    ) +
    ggplot2::ggtitle(title) + ggplot2::ylab("Scaled Euclidean distance") +
    ggplot2::scale_y_reverse(limits = rev(ylim)) +
    ggplot2::theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    ) +
    ggplot2::scale_x_discrete(labels = x_lables)
  return(list(plot_object))
}

#' Plot multiple ggplot objects to one PDF
#'
#' @param plot_objects A list of ggplot objects.
#'
#' @param file File path of the PDF to be created.
#'
#' @param nrow Number of rows of plot objects in one page.
#'
#' @param ncol Number of columns of plot objects in one page.
#'
#' @param width Width of the plot size.
#'
#' @param height Height of the plot size.
#'
#' @param units Units("in", "cm", "mm", or "px") of the plot size.
#'
ggplot_to_pdf <-
  function(plot_objects, file, nrow, ncol, width, height, units) {
    ggplot2::ggsave(
      filename = file,
      plot = gridExtra::marrangeGrob(plot_objects, nrow = nrow, ncol = ncol),
      width = width, height = height, units = units
    )
  }

#' Get signature activity information from an exposure matrix
#'
#' @param exposure Exposures as a numerical matrix (or data.frame) with
#'   signatures in rows and samples in columns. Rownames are taken as the
#'   signature names and column names are taken as the sample IDs.
#'
#' @return A data frame showing the activity information for each signature in
#'   \code{exposure}.
#'
get_sig_activity <- function(exposure) {
  get_sig_info <- function(sig_exposures) {
    total_tumors <- length(sig_exposures)
    active_tumors <- sum(sig_exposures > 0)
    sig_prop <- active_tumors / total_tumors
    return(list(
      active_tumors = active_tumors,
      sig_prop = sig_prop
    ))
  }
  
  sig_info <- apply(exposure, MARGIN = 1, FUN = get_sig_info)
  active_tumors <- sapply(sig_info, FUN = "[[", 1)
  sig_prop <- sapply(sig_info, FUN = "[[", 2)
  
  sig_activity_df <- data.frame(
    sig_id = names(active_tumors),
    active_tumors = active_tumors,
    sig_prop = sig_prop
  )
  
  sig_activity_df <-
    dplyr::arrange(sig_activity_df, dplyr::desc(active_tumors))
  
  return(sig_activity_df)
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
multiple_boxplots <- function(distance_df, data_type, ylim) {
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

