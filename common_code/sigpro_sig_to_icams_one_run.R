# input_path: path to store SigPro-formatted signature.
# e.g.
# SBS_set1/raw_results/SigProfilerExtractor.results/
# Realistic/seed.145879/SBS96/Suggested_Solution/
# SBS96_De-Novo_Solution/Signatures/SBS96_De-Novo_Signatures.txt
#
# output_path (optional): The directory to save ICAMS-formatted 
# CSV signature file. Usually refers to the grand-parent directory 
# of input_path with basename "seed.<...>", e.g.
# SBS_set1/raw_results/SigProfilerExtractor.results/
# Realistic/seed.145879/
# It can be other directory, however.
#
# cat_type: "SBS96" or "ID83"
sigpro_sig_to_icams_one_run <- function(
    input_path, 
    output_path = paste0(dirname(input_path), 
                         "/../../../../"),
    cat_type) {
  #browser() # debug
  stopifnot(cat_type %in% c("SBS96", "ID83"))
  sig_catalog_sp <- utils::read.table(
    input_path,
    sep = "\t",
    as.is = TRUE,
    header = TRUE)
  if (cat_type == "SBS96") {
    sig_catalog <- ICAMS:::MakeSBS96CatalogFromSigPro(sig_catalog_sp)
  } else if (cat_type == "ID83") {
    sig_catalog <- ICAMS:::MakeID83CatalogFromSigPro(sig_catalog_sp)
  }
  sig_catalog <- ICAMS::as.catalog(sig_catalog,
                                   catalog.type = "counts.signature")
  if(!is.null(output_path) && !is.na(output_path)) {
    ICAMS::WriteCatalog(sig_catalog, 
                        file.path(output_path, "extracted.signatures.csv"))
  }
  invisible(return(sig_catalog))
}