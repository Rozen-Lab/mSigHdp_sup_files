# Code from Steve Rozen, 2021 01 24
# Analyses for Liu Mo's paper, "mSigHdp: hierarchical Dirichlet processes
# in mutational signature extraction, Liu et al."
#
# This code assembles run-level results that are in individual .csv files
# into one large tibble - all_results.csv,
# and then does two levels of summary:
#
# stats_by_approach_and_data_set.csv - at the level of approach and data set, 
# stats_by_approach.csv - at the level of approach

# Please run this script from the top directory
if (basename(getwd()) != "Liu_et_al_Sup_Files") {
  stop("Please run from top level directory, Liu_et_al_Sup_Files")
}

# 0. Install and load dependencies --------------------------------------------

library(data.table)
library(dplyr)
library(readr)
library(tibble)



# 1. Combine all extraction measures ------------------------------------------

myread <- function(ff) as_tibble(fread(ff)[ , -1])

# These are the files to be logically "cbind"ed
gtSigs <- paste0("ID",c(1:6,8,9,11,13,14))
nSigs <- length(gtSigs)


home <-
  "other_analyses/indel_beta_selection/summary/top_level_summary/"

fff <- c("averCosSim.csv",
         "PPV.csv",
         "TPR.csv",
         "NumSigsExtracted.csv",
         "falseNeg.csv",
         "falsePos.csv")

dd <- lapply(paste0(home, "/", fff), myread)

# Clean up all csv files after reading
unlink(paste0(home, "/*.csv"))


# Check to make sure the rows in all the files were
# for the same approach, data-set, and run
com <- dd[[1]]
## Check the 1st, 3rd and 4th columns in the files.
for (ii in c(1, 3:4)) {
  cat(
    all(
      unlist(
        lapply(dd, function(xx, col) all(com[ ,  col] == xx[ , col]), ii))),
    "\n")
}

# Make one big tibble (cc)
cc <- dd[[1]]
for (ii in 2:length(dd)) {
  cc2 <- full_join(cc, dd[[ii]])
  cc <- cc2
}

# Make the column names more user-friendly
ccc <- cc[ , c(3, 4, 1, 2, 5:ncol(cc))]
colnames(ccc) <- c("Approach",
                   "Noise_level",
                   "Run",
                   "aver_Sim_TP_only",
                   "PPV",
                   "TPR",
                   "N_sigs",
                   "FN",
                   "FP")
# Add Composite Measure after column "Run".
ccc <- ccc %>% 
  mutate(Composite = PPV + TPR + aver_Sim_TP_only,
         .after = Run)

readr::write_csv(ccc, paste0(home, "/../all_results.csv"))




# 2. Summarize extraction measures by software packages and data set ----------

by.approach.and.data.set <- ccc %>%
  group_by(Approach, Noise_level) %>%
  summarise_at(.vars = colnames(ccc)[4:ncol(ccc)], .funs = c(mean, sd))

# Name-fixing function for summary tibbles with two-statistics (mean, sd)
fix.names <- function(nn) {
  vv <- sub("(.*)_fn1", "mean(\\1)", x = nn, perl = TRUE)
  sub("(.*)_fn2", "sd(\\1)", x = vv, perl = TRUE)
}

# Name-fixing function for summary tibbles with three-statistics (mean, median sd)
fix.names.3stats <- function(nn) {
  vv <- sub("(.*)_fn1", "mean(\\1)", x = nn, perl = TRUE)
  ww <- sub("(.*)_fn2", "median(\\1)", x = vv, perl = TRUE)
  sub("(.*)_fn3", "sd(\\1)", x = ww, perl = TRUE)
}


colnames(by.approach.and.data.set) <- fix.names(colnames(by.approach.and.data.set))
#num.measures <- 7 + 2 * nSigs
num.measures <- 7
perm <- c(1:2, unlist(lapply(3:(2+num.measures), function(x) c(x, x+num.measures))))
colnames(by.approach.and.data.set)[perm] # check
by.approach.and.data.set <- by.approach.and.data.set[ , perm]

readr::write_csv(by.approach.and.data.set, 
                 paste0(home, "stats_by_approach_and_data_set.csv"))



# 3. Summarize by approach ----------------------------------------------------
#
# This one is for checking supplementary tables, which were assembled
# separately.  Summarize by approach, and re-order / rename the
# columns
by.approach <- ccc %>%
  group_by(Approach) %>%
  summarise_at(.vars = colnames(ccc)[4:ncol(ccc)], .funs = c(mean, median, sd))
colnames(by.approach) <- fix.names.3stats(colnames(by.approach))


perm3 <- c(1, unlist(lapply(2:(1+num.measures), function(x) c(x, x+num.measures, x+2*num.measures))))
colnames(by.approach)[perm3] # check
by.approach <- by.approach[ , perm3]

readr::write_csv(by.approach, paste0(home, "stats_by_approach.csv"))

# 4. Wilcoxon rank sum test between approaches --------------------------------
#
# This test is to evaluate whether the composite measure between approaches
# are significantly different.

approaches <- ccc %>% select("Approach") %>% unique() %>% unlist() %>% unname()
index <- length(approaches)

pairwise.signif <- data.frame(
  approach1 = character(0),
  approach2 = character(0),
  p.value = numeric(0))


for(ii in seq(1,index-1)){
  
  for (jj in seq(ii+1,index)){
    
    comp <- ccc %>%
      select(Approach,Composite)
    
    compII <- comp %>%
      filter(Approach %in% approaches[ii]) %>%
      select(Composite) %>%
      unlist() %>% unname()
    
    compJJ <- comp %>%
      filter(Approach %in% approaches[jj]) %>%
      select(Composite) %>%
      unlist() %>% unname()
    
    res <- stats::wilcox.test(compII,compJJ)
    
    current <- data.frame(
      approach1 = approaches[ii],
      approach2 = approaches[jj],
      p.value = res$p.value)
    
    pairwise.signif <- rbind(pairwise.signif,current)
    
  }
  
}

readr::write_csv(pairwise.signif, paste0(home, "Pairwise.Wilcoxon.Rank.Sum.csv"))

