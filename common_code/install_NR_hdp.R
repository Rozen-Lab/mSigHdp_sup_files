# Install and load the specific package versions of hdpx and mSigHdp 
# needed to test Nicola Roberts's algorithms.

# The specific package versions must be installed before running
# any data generate that uses common_code/generic_run_NR_hdp.R.

# WARNING: sourcing this file will remove other versions of hdpx and mSigHdp.

# Please contact Steve Rozen, steverozen@pm.me, if you have difficulties.

if (!file.exists("common_code/set_top_dir.R")) {
  stop("Please run from top level directory")
}

library(remotes)

hdpx.version <- "0.1.5.0099"
if (system.file(package = "hdpx") != "") {
  # Re-install if hdpx fails to be loaded, 
  # or its version does not match hdpx.version
  if (!requireNamespace("hdpx", quietly = TRUE) ||
      packageVersion("hdpx") != hdpx.version) {
    remove.packages("hdpx")
    remotes::install_github("steverozen/hdpx", ref = "NR-version-plus-fixes")
  }
} else {
  # Install directly if not been installed
  remotes::install_github("steverozen/hdpx", ref = "NR-version-plus-fixes")
}
message("hdpx version ", packageVersion("hdpx"))
stopifnot(packageVersion("hdpx") == hdpx.version)

mSigHdp.version <- "0.0.0.9019"
# Re-install if mSigHdp fails to be loaded, 
# or its version does not match mSigHdp.version
if (system.file(package = "mSigHdp") != "") {
  if (!requireNamespace("mSigHdp", quietly = TRUE) || 
      packageVersion("mSigHdp") != mSigHdp.version) {
    remove.packages("mSigHdp")
    remotes::install_github(repo = "steverozen/mSigHdp", 
                            ref = "for-NR-version-plus-fixes")
  }
} else {
  remotes::install_github(repo = "steverozen/mSigHdp", 
                          ref = "for-NR-version-plus-fixes")
}
message("mSigHdp version ", packageVersion("mSigHdp"))
stopifnot(packageVersion("mSigHdp") == mSigHdp.version)
