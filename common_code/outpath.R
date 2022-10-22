# Function to generate a path under folder "output_for_paper"
outpath <- function(filename) {
  if(!dir.exists("output_for_paper")) dir.create("output_for_paper")
  file.path("output_for_paper", filename)
}