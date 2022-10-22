# Please run this script from the top directory
if (basename(getwd()) != "Liu_et_al_Sup_Files") {
  stop("Please run from top level directory, Liu_et_al_Sup_Files")
}

# Install and load dependencies -----------------------------------------------

require(ICAMS)
require(SynSigEval)



# Specify global options ------------------------------------------------------
options(stringsAsFactors = F)



# Specify global variables ----------------------------------------------------

topLevelFolder4Data <- "./indel/input/"
topLevelFolder4Run <-
  "./other_analyses/indel_beta_selection/raw_results"
folder4ToolwiseSummary <- 
  "./other_analyses/indel_beta_selection/summary/toolwise_summary"
folder4CombinedSummary <- 
  "./other_analyses/indel_beta_selection/summary/top_level_summary"

# Specify name of ground-truth spectra data set
datasetNames <- "Realistic"

# Specify names of computational approaches
# R-based tools which can do both extraction and attribution.
#
# SigProfilerExtractor is included, because its extracted signature file is
# already transformed into ICAMS-csv format when running 
# "code/2_ID_eval/2_SP_tsv_to_csv.R"
RBasedExtrAttrToolNames <- c(paste0("mSigHdp.beta",c(40, 50, 60, 75)),
                             "SigProfilerExtractor")

# Specify seeds used in analysis.
# Specify 5 seeds used in software running
seedsInUse <- c(145879, 200437, 310111, 528401, 1076753)



# Summarize on individual runs ------------------------------------------------

for (datasetName in datasetNames) {
  for (seedInUse in seedsInUse) {
    ## Summarize R-based Extraction and attribution tools.
    for (extrAttrToolName in RBasedExtrAttrToolNames) {
      SynSigEval::SummarizeSigOneExtrAttrSubdir(
        run.dir = paste0(topLevelFolder4Run, "/",
                         extrAttrToolName, ".results/", 
                         datasetName, "/seed.",seedInUse,"/"),
        ground.truth.exposure.dir = paste0(topLevelFolder4Data,"/",
                                           datasetName,"/"),
        summarize.exp = F,
        overwrite = T
      )
    }
  }
}




# Summary of runs on each dataset with each software --------------------------

for (datasetName in datasetNames) {
  ## For each dataset, summarize 20 runs
  ## using different seeds by EMu
  for (extrAttrToolName in RBasedExtrAttrToolNames) {
    SynSigEval::SummarizeMultiRuns(
      datasetName = datasetName,
      toolName = extrAttrToolName,
      resultPath = paste0(topLevelFolder4Run, "/",
                          extrAttrToolName, ".results/",
                          datasetName, "/"),
      run.names = paste0("seed.",seedsInUse)
    )
  }
}



# Summarize results of multiple data sets by each tool ------------------------

datasetGroup <- c("noisy.NB.size.10")
names(datasetGroup) <- datasetNames

toolsToEval <- RBasedExtrAttrToolNames

for(toolName in toolsToEval){
  SummarizeOneToolMultiDatasets(
    datasetNames = datasetNames,
    datasetGroup = datasetGroup,
    datasetGroupName = "Noise level",
    datasetSubGroup = NULL,
    datasetSubGroupName = NULL,
    toolName = toolName,
    toolPath = paste0(topLevelFolder4Run,"/",toolName,".results/"),
    out.dir = paste0(folder4ToolwiseSummary,"/",toolName,"/"),
    display.datasetName = FALSE,
    overwrite = T
  )
}



# Summarize results of multi tools on multi data sets -------------------------

FinalExtrAttr <- SummarizeMultiToolsMultiDatasets(
  toolSummaryPaths = paste0(folder4ToolwiseSummary,"/",toolsToEval,"/"),
  out.dir = folder4CombinedSummary,
  display.datasetName = FALSE,
  sort.by.composite.extraction.measure = "descending",
  overwrite = T
)


