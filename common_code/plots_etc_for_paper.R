if (!file.exists("common_code/set_top_dir.R")) {
  stop("Please run from top level directory")
}

# Loading required packages ---------------------------------------------------
library(data.table)
library(tibble)
library(magrittr)
library(dplyr)
library(beeswarm)
library(MASS)
library(tidyr)
library(openxlsx) # https://cran.r-project.org/web/packages/openxlsx/openxlsx.pdf


# Utility functions -----------------------------------------------------------
source("common_code/outpath.R")

set1_pch <- 16
set2_pch <- 17

split_by_approach_and_pull <- function(vv, approaches.to.use) {
  xxs     <- split(vv, vv$Approach)

  # This makes sure that the elements of xx2 are in the same order as 
  # the elements of approaches.to.use, which should be in the order 
  # we want for columns of the eventual beeswarm graph.
  xxs2    <- xxs[approaches.to.use] 
  
  my.pull <- function(my.approach, colname) { # only use my.pull for debugging
                                              # anly apply it to approaches.to.use
    zz <- xxs[[my.approach]]
    if (length(zz) == 0) {
      stop("No results for my.approach = ", my.approach, " colname = ", colname)
      
    }
    pull(zz, colname)}
  comp <- lapply(xxs2, pull, "Composite")
  tpr  <- lapply(xxs2, pull, "TPR")
  ppv  <- lapply(xxs2, pull, "PPV")
  sim  <- lapply(xxs2, pull, "aver_Sim_TP_only")
  return(list(split = xxs2, comp = comp, tpr = tpr, ppv = ppv, sim =sim))
}


# Definition of plotting functions --------------------------------------------
four_beeswarms <- function(ww, main, col, pch, filename,
                           mfrow = c(3,1),
                           legend.fn = NULL,
                           mar = c(9, 12, 4, 12) + 0.1,
                           override.labels = NULL) {

  grDevices::cairo_pdf(
    filename = outpath(filename),
    height = 9, 
    onefile = TRUE)
  par(mfrow = mfrow, mar = mar)

  beeswarm(x = ww$comp, las = 2, ylab = "Composite Measure", 
           pwpch = pch, pwcol = col, labels = "")
  if (!is.null(legend.fn)) {
    legend.fn()
  }
  
  beeswarm(x = ww$ppv, las = 2,  ylab = "PPV", 
           pwpch = pch, pwcol = col, labels = "")
  
  beeswarm(x = ww$tpr, las = 2, ylab = "TPR", 
           pwpch = pch, pwcol = col, labels = "")
  
  if (!is.null(override.labels)) {
  beeswarm(x = ww$sim, las = 2, ylab = "Cosine similarity", 
           pwpch = pch, pwcol = col, labels = override.labels)
  } else {
    beeswarm(x = ww$sim, las = 2, ylab = "Cosine similarity", 
             pwpch = pch, pwcol = col)
  }

  grDevices::dev.off()
  
}


generic_4_beeswarm_fig <-
  function(tt, 
           approaches.to.use, # character vector of names of approaches
           sbs.or.indel, 
           file.name.prefix,
           legend.fn = NULL,
           col       = "black",
           mfrow = c(3, 1),
           mar = c(8, 15, 4, 14) + 0.1,
           override.labels = NULL) {

  set1 <- paste0(sbs.or.indel, "_set1")
  set2 <- paste0(sbs.or.indel, "_set2")
  t1 <- filter(tt, Noise_level == "Realistic")
  xx <- filter(t1, Data_set %in% c(set1, set2))

  ww <- split_by_approach_and_pull(xx, approaches.to.use)
  xx.data.set <- unlist(lapply(ww$split, pull, "Data_set"))

  col <- ifelse(xx.data.set == set1, "black", "black")
  pch <- ifelse(xx.data.set == set1, set1_pch, set2_pch)
  
  four_beeswarms(ww,
                 col             = col,
                 pch             = pch,
                 filename        = paste0(file.name.prefix, sbs.or.indel, ".pdf"),
                 legend.fn       = legend.fn,
                 mfrow           = mfrow,
                 mar             = mar,
                 override.labels = override.labels)
}


downsample.approaches <- c("mSigHdp", 
                           "mSigHdp_ds_10k", 
                           "mSigHdp_ds_5k",
                           "mSigHdp_ds_3k", 
                           "mSigHdp_ds_1k")

downsample.override <- c(
  "None", "10,000", "5,000", "3,000", "1,000"
)


downsample_indel_fig <- function(tt) {
  approaches <- downsample.approaches
  generic_4_beeswarm_fig(tt                = tt, 
                         approaches.to.use = approaches, 
                         sbs.or.indel      = "indel", 
                         file.name.prefix  = "sup_fig_s12_downsampling_",
                         mfrow = c(3, 1),
                         mar = c(8, 14, 4, 14) + 0.1,
                         legend.fn = function() { 
                           set1_set2_legend("indel")},
                         override.labels = downsample.override)
  
}


downsample_SBS_fig <- function(tt) {
  approaches <- downsample.approaches
  generic_4_beeswarm_fig(tt, approaches, "SBS", "sup_fig_s7_downsampling_",
                         mfrow = c(3, 1),
                         mar = c(8, 14, 4, 14) + 0.1,
                         legend.fn = function() { 
                           set1_set2_legend("SBS", "topleft")
                         },
                         override.labels = downsample.override)
  
}


noise_level_fig <- function(tt, indel.or.sbs, approach) {

  # filter by approach and by indel_set1
  # split by approach, color by noise level
  data.set = paste0(indel.or.sbs, "_set1")
  tt1 <- dplyr::filter(tt, Data_set == data.set)
  ww <- split_by_approach_and_pull(tt1, approach)
  xx.noise.level <- unlist(lapply(ww$split, pull, "Noise_level"))
  
  # Color by noise
  col <- ifelse(xx.noise.level == "Realistic", "red",
                ifelse(xx.noise.level == "Moderate", "violet", "blue"))
  
  pch <- ifelse(xx.noise.level == "Realistic", set1_pch,
                ifelse(xx.noise.level == "Moderate", set1_pch, set1_pch))
  
  if (indel.or.sbs == "indel") {
    filename <- "sup_fig_s15_indel_noise.pdf"
  } else if (indel.or.sbs == "SBS") {
    filename <- "sup_fig_s11_SBS_noise.pdf"
  } else {
    stop("indel.or.sbs expects to be 'indel' or 'SBS'")
  }
  
  four_beeswarms(
    ww,
    col      = col,
    pch      = pch,
    filename = filename,
    legend.fn = function() {
      legend(x = "bottomleft",
             title  = "Resampling noise",
             legend = c("None", "Moderate", "Realistc"),
             col    = c("blue",  "violet", "red"),
             pch    = 18:16,
             bty    = "n")
    })

}


SBS35_detect <- function(tt) {

  sigpro.pch  <- "S"
  msighdp.pch <- "M"
  ds1.col     <- "cornflowerblue"
  ds2.col     <- "coral"
  
  spike.in.counts <- c(5L, 10L, 20L, 30L, 50L, 100L)
  data.sets.1066 <- paste0("sens_SBS35_", spike.in.counts, "_1066")
  data.sets.728 <- paste0("sens_SBS35_", spike.in.counts, "_728")
  data.sets <- c(data.sets.1066, data.sets.728)
  data.set.2.count <- rep(spike.in.counts, 2)
  names(data.set.2.count) <- data.sets
  tt1 <- dplyr::filter(tt, Data_set %in% data.sets)
  tt2 <- dplyr::filter(tt1, Noise_level == "Realistic")
  tt3 <- dplyr::filter(tt2, 
                       Approach %in% c("mSigHdp_ds_3k",
                                       "SigProfilerExtractor"))
  FNs <- dplyr::pull(tt3, FN.sigs)
  SBS35.found <- unlist(lapply(FNs, function(zz) !any(grepl("SBS35", zz, fixed = TRUE))))
  tt4 <- dplyr::mutate(tt3, SBS35.found = SBS35.found)
  tt5 <- tt4[ , c("Data_set", "Approach", "Run", "SBS35.found", "FN.sigs")]
  dplyr::group_by(tt5, Data_set, Approach) %>% 
    dplyr::summarise(num.found = sum(SBS35.found), .groups = "drop") -> tt6

  sbs35_detect <- mutate(tt6, spike.in.count = data.set.2.count[Data_set])
  
  # rrx <- robust::lmRob(formula = num.found ~ spike.in.count + as.factor(Approach),
  #                      data = sbs35_detect)
  
  data.table::fwrite(sbs35_detect, file = outpath("sup_fig_s10_sbs35_detect_data.csv"))
  
  df <- nrow(sbs35_detect) - 3
  sens_stats_r <- 
    MASS::rlm(formula = num.found ~ spike.in.count + Approach,
              data = sbs35_detect)
  sens_stats_rs <- summary(sens_stats_r)
  sens_coef_r <- sens_stats_rs$coefficients
  sens_stats_rp <- 2*pt(-abs(sens_stats_rs$coefficients[, 3]), df = df)
  sbs35_detect_coef <- cbind(sens_coef_r, p = sens_stats_rp)
  sbs35_detect_coef <- data.frame(sbs35_detect_coef)
  sbs35_detect_coef <- cbind(variable = rownames(sbs35_detect_coef), sbs35_detect_coef)
  openxlsx::write.xlsx(sbs35_detect_coef,
                       file = outpath("sup_fig_s10_sbs35_detect_coef.xlsx"))
    
  sens_stats_l <-
    lm (formula = num.found ~ spike.in.count + Approach, 
        data = sbs35_detect)
  sens_stats_ls <- summary(sens_stats_l)
  sens_stats_lp <- 2*pt(-abs(sens_stats_ls$coefficients[ , 3]), df = df)
  
  capture.output(df, sens_stats_rs, sens_stats_rp,
                 sens_stats_ls, sens_stats_lp,
                 file = outpath("sup_fig_s10_sbs35_detect_stats.txt"))
  
  save(sens_stats_r, sens_stats_l, sbs35_detect, 
       file = outpath("sup_fig_s10_sbs35_detect.Rdata"))
  
  to.plot <- split(sbs35_detect, sbs35_detect$spike.in.count)
  to.plot <- to.plot[as.character(spike.in.counts)]
  to.plot2 <- lapply(to.plot, pull, num.found)
  to.plot2.app <- unlist(lapply(to.plot, pull, Approach))
  to.plot2.dset <- unlist(lapply(to.plot, pull, Data_set))

  pch <- ifelse(to.plot2.app == "mSigHdp_ds_3k", msighdp.pch, sigpro.pch)

  grDevices::cairo_pdf(
    filename = outpath("sup_fig_s10_sbs35_detect.pdf"),
    height = 4, 
    onefile = TRUE)
  par(mar = c(5.1, 5.1, 4.8, 2.1), xpd = TRUE)

  beeswarm(x = to.plot2, las = 2, 
           ylab = "Number of runs in each data\nset in which SBS35 was detected", 
           xlab = "Number of synthetic spectra containing SBS35",
           pwpch = pch, spacing = 1.6)
  legend.info <- legend(
    x         = 0,5,
    y         = 7.4,
    title     = "Legend",
    title.adj = 0,
    legend    = 
      c("mSigHdp run on one of the 2 datasets at a given x-axis position",
        "SigProfilerExtractor run on one of the 2 datasets at a given x-axis position"),
    pch       = c(msighdp.pch, sigpro.pch),
    bty       = "n",
    col       = "black",
    border   = "white",
    fill     = NULL,
    lty      = 0
  )
  dev.off()
  
  invisible(to.plot)
}


main_text_cpu <- function(sbs.or.indel, approaches.to.use) {
  uu <- data.table::fread(outpath("supplementary_table_s5.csv"))
  data.sets <- paste0(sbs.or.indel, "_set", c(1, 2))
  filter(uu, 
         Data_set %in% data.sets &
           Approach %in% approaches.to.use) %>%
    mutate(CPU.days = cpu_time / (60 * 60 * 24)) %>%
    mutate(Noise_level = NULL, cpu.time = NULL) -> uu3
  # browser()
  to.plot <- split(uu3, uu3$Approach)
  to.use <- which(approaches.to.use %in% names(to.plot))
  approaches.to.use <- approaches.to.use[to.use]
  to.plot <- to.plot[approaches.to.use]

  dplyr::group_by(uu3, Data_set, Approach) %>%
    dplyr::summarise(mean_CPU_days = mean(CPU.days)) -> return.value
  
  to.plot2 <- lapply(to.plot, pull, CPU.days)
  
  xx.data.set <- unlist(lapply(to.plot, pull, "Data_set"))
  
  col <- ifelse(xx.data.set == data.sets[1], "black", "black")
  pch <- ifelse(xx.data.set == data.sets[1], set1_pch, set2_pch)
  
  beeswarm(x      = to.plot2, 
           las    = 2, 
           ylab   = "CPU days",
           pwpch  = pch,
           pwcol = col)
  
  legend(x = "topright",
         legend = paste0(sbs.or.indel, "_set", 1:2),
         col    = "black", # c("red",  "blue"),
         pch    = c(16,     17),
         bty    = "n")
  
  invisible(return.value)
}

set1_set2_legend <- function(sbs.or.indel, pos = "bottomleft") {
  legend(x = pos,
         legend = paste0(sbs.or.indel, "_set", 1:2),
         col    = c("black",  "black"),
         pch    = c(16,     17),
         bty    = "n")
}


main_text_table <- function(tt, approaches.to.use, sbs.or.indel) {
  set1 <- paste0(sbs.or.indel, "_set1")
  set2 <- paste0(sbs.or.indel, "_set2")

  t1 <- filter(
    tt, 
    Data_set %in% c(set1, set2), 
    Noise_level == "Realistic",
    Approach %in% approaches.to.use) -> t2
  
  if (sbs.or.indel == "SBS") {
    best <- "mSigHdp_ds_3k"
  } else {
    best <- "mSigHdp"
  }
  wt <- wilcox.test(Composite ~ Approach, 
                    data = t2, 
                    subset = Approach %in% c(best, "SigProfilerExtractor"))

  if (sbs.or.indel == "SBS") {
    wt2 <- wilcox.test(Composite ~ Approach, 
                      data = t2, 
                      subset = Approach %in% c(best, "mSigHdp"))
  }

  dplyr::group_by(t2, Data_set, Approach) %>%
    dplyr::summarise(mean_comp = mean(Composite),
                     sd_comp   = sd(Composite),
                     mean_PPV  = mean(PPV),
                     mean_TPR  = mean(TPR),
                     mean_cos  = mean(aver_Sim_TP_only)
    ) %>%
    arrange(desc(mean_comp), .by_group = TRUE) -> t3
  
  
  t3 %>% filter(Data_set == set1) %>% nrow -> num.set1
  # browser()
  colnames(t3) <- c("Data\nset", "Approach",
                    "Mean", "SD", # For composite measure
                    "Mean\nPPV", "Mean\nTPR",
                    "Mean cosine\nsimilarity")
  
  wb <- openxlsx::createWorkbook()
  
  heading.style <- 
    createStyle(halign = "center", textDecoration = "bold", 
                wrapText = TRUE, valign = "center",
                borderStyle = "medium", border = "bottom")
  
  num.style <- 
    createStyle(numFmt = "0.00", halign = "center")
  
  left.style <- 
    createStyle(valign = "center", halign = "center", 
                textDecoration = "bold")
  
  top.border.style <- createStyle(border = "top", borderStyle = "medium")
  
  addWorksheet(wb, sbs.or.indel)
  
  startrow <- 1
  writeData(wb, 1, "Composite Measure", startCol = 3, startRow = startrow)
  mergeCells(wb, 1, cols = 3:4, rows = startrow) # Merge "Composite measure"
  
  datarow1 <- startrow + 2
  writeData(wb, 1, t3, startRow = datarow1 - 1, startCol = 1)
  set2.datarow1 <- datarow1 + num.set1
  
  mergeCells(wb, 1, cols = 1, rows = datarow1:(set2.datarow1 - 1))

  last.datarow <- datarow1 + nrow(t3) - 1
  
  mergeCells(wb, 1, cols = 1, rows = set2.datarow1:last.datarow)
   
  addStyle(wb, 1, num.style, cols = 3:7,  
           rows = datarow1:last.datarow,
           gridExpand = TRUE)
  
  addStyle(wb, 1, heading.style, 
           cols = 1:7, rows = startrow + 1, gridExpand = TRUE)
  addStyle(wb, 1, 
           heading.style, cols = 3:4, rows = startrow, gridExpand = TRUE)

  addStyle(wb, 1, left.style, cols = 1, rows = datarow1:last.datarow)
  
  addStyle(wb, 1, top.border.style, cols = 1:7, rows = last.datarow + 1)
  
  writeData(wb, 1, startCol = 1, startRow = last.datarow + 2,
            paste0("Wilcoxon rank-sum test ", best, " vs ", "SigProfilerExtractor p = ",
            format(wt$p.value, scientific = TRUE, digits = 4)))
  if (sbs.or.indel == "SBS") {
    writeData(wb, 1, startCol = 1, startRow = last.datarow + 3,
              paste0("Wilcoxon rank-sum test ", best, " vs ", "mSigHdp p = ",
                     format(wt2$p.value, scientific = TRUE, digits = 4)))
  }

  if (sbs.or.indel == "SBS") {
    filename <- "table_1_SBS.table.xlsx"
  } else if (sbs.or.indel == "indel") {
    filename <- "table_2_indel.table.xlsx"
  } else {
    stop("sbs.or.indel expects to be 'SBS' or 'indel'")
  }
  
  openxlsx::saveWorkbook(wb, 
                         outpath(filename), 
                         overwrite = TRUE)
  
  t3
}


cpu_fig_and_table <- function() {
  
  grDevices::cairo_pdf(
    filename = outpath("fig_7_CPU_time.pdf"),
    height   = 11,
    width    = 7, 
    onefile = TRUE)
  par(mfrow = c(2, 1), mar = c(9, 8, 4, 8) + 0.1)
  
  sbs.cpu   <- main_text_cpu("SBS",   c("mSigHdp_ds_3k",
                                        "mSigHdp",
                                        "SigProfilerExtractor",
                                        "signeR",
                                        "SignatureAnalyzer"))
  
  indel.cpu <- main_text_cpu("indel", c("mSigHdp",
                                        "SigProfilerExtractor",
                                        "signeR",
                                        "SignatureAnalyzer"))
  
  dev.off()
  
  cpu.summary <- rbind(sbs.cpu, indel.cpu)
  tidyr::pivot_wider(cpu.summary, 
                     names_from = Data_set, 
                     values_from = mean_CPU_days) %>%
    mutate(`Average SBS` = (SBS_set1 + SBS_set2) / 2, 
           .keep = "all", 
           .before = "indel_set1") %>%
    mutate(`Average indel` = (indel_set1 + indel_set2) / 2,
           .keep = "all") -> tmp.table
    tmp.table <- tmp.table[c(2, 1, 5, 4, 3)  , ]
    openxlsx::write.xlsx(tmp.table, outpath("table_3_cpu_time.xlsx"))
}


# Create either the SBS or the indel part of sup table s6
sup_table_s6 <- function(tt, sbs.or.indel) {
  dplyr::filter(tt,
                Approach %in% c("SigProfilerExtractor",
                                "SP_Kplus2",
                                "SP_KmSigHdp") &
                  Data_set == paste0(sbs.or.indel, "_set1") &
                  Noise_level == "Realistic") %>%
    dplyr::transmute(Approach  = Approach,
                     seed      = sub("seed\\.", "", Run),
                     N_Sigs    = N_Sigs,
                     FN        = FN,
                     FP        = FP,
                     Composite = Composite) -> xx
  
  # Make 3 tables
  dplyr::filter(xx, Approach == "SigProfilerExtractor") %>% dplyr::mutate(Approach = NULL) -> sp
  dplyr::filter(xx, Approach == "SP_Kplus2") %>% dplyr::mutate(Approach = NULL, seed = NULL) -> plus2
  dplyr::filter(xx, Approach == "SP_KmSigHdp") %>% dplyr::mutate(Approach = NULL, seed = NULL) -> ms
  
  add_increases <- function(uu) {
    dplyr::mutate(uu, TP_increase  = sp$FN - FN, .after = FN) %>%
      dplyr::mutate(FP_increase = FP - sp$FP, .after = FP)
  }
  
  plus2 <- add_increases(plus2)
  ms    <- add_increases(ms)
  
  ms <- dplyr::mutate(ms, N_Sigs = NULL)
  colnames(plus2) <- paste("plus2", colnames(plus2), sep = "_")  
  colnames(ms) <- paste("msigK", colnames(ms), sep = "_")
  retval <- cbind(sp, plus2, ms)[c(2:5,1), ]
  
  plus2_wt <- wilcox.test(retval$Composite, retval$plus2_Composite)
  
  msigK_wt <- wilcox.test(retval$Composite, retval$msigK_Composite)
  
  avg <- cbind(seed = "Mean", as_tibble(t(colMeans(retval[ , 2:length(retval)]))))
  wtp <- t(c("P-value", rep("", 9), plus2_wt$p.value, rep("", 4), msigK_wt$p.value))
  colnames(wtp) <- colnames(retval)
  retval <- rbind(retval, avg, wtp)
  
  pr <- outpath(paste0("sup_table_s6_", sbs.or.indel))
  data.table::fwrite(retval, paste0(pr, ".csv"))
  invisible(retval)
}


all_figs_and_tables_this_file <- function(tt) {
  # tt should be the output of summarize_all_level1_dirs in file summarize_level1_dirs.R
  
  main.text.SBS.approaches <- 
    c("mSigHdp_ds_3k",
      "mSigHdp",
      "SigProfilerExtractor",
      "NR_hdp_gb_20",
      "NR_hdp_gb_1",
      "signeR",
      "SignatureAnalyzer")

  generic_4_beeswarm_fig(
    tt = tt, 
    approaches.to.use = main.text.SBS.approaches,
    sbs.or.indel      = "SBS", 
    file.name.prefix  = "fig_5_main_text_",
    mfrow             = c(3, 1),
    mar               = c(8, 14, 4, 14) + 0.1,
    legend.fn         = function() { set1_set2_legend("SBS")})
  
  main_text_table(
    tt                = tt,
    approaches.to.use = main.text.SBS.approaches,
    sbs.or.indel      = "SBS"
  )
  
  main.text.indel.approaches <- # Important, this is the order in the 
                                # main text figure; make sure it is 
                                # sorted by the average composite measure
                                # in indel_set1 and indel_set2
    c("mSigHdp",
      "NR_hdp_gb_50", 
      "NR_hdp_gb_1",
      "SigProfilerExtractor",
      "SignatureAnalyzer",
      "signeR")
  
  main_text_table(
    tt                = tt,
    approaches.to.use = main.text.indel.approaches,
    sbs.or.indel      = "indel"
  )
  
  generic_4_beeswarm_fig(
    tt = tt, 
    approaches.to.use = main.text.SBS.approaches,
    sbs.or.indel      = "SBS", 
    file.name.prefix  = "fig_5_main_text_",
    mfrow             = c(3, 1),
    mar               = c(8, 14, 4, 14) + 0.1,
    legend.fn         = function() { set1_set2_legend("SBS")})
  
  generic_4_beeswarm_fig(
    tt = tt,
    approaches.to.use = main.text.indel.approaches, 
    sbs.or.indel      = "indel",
    file.name.prefix  = "fig_6_main_text_",
    mfrow = c(3, 1),
    mar = c(8, 14, 4, 14) + 0.1,
    legend.fn = function() { set1_set2_legend("indel")})
  
  noise_level_fig(
    tt, 
    "indel",
    approach = c("mSigHdp",
                 "SigProfilerExtractor",
                 "SignatureAnalyzer",
                 "signeR")) # Order of SA and signeR are reversed between indel and SBS
  
  
  noise_level_fig(tt, "SBS",approach = c("mSigHdp_ds_3k",
                                         "mSigHdp",
                                         "SigProfilerExtractor",
                                         "signeR",
                                         "SignatureAnalyzer"))
  downsample_indel_fig(tt)
  downsample_SBS_fig(tt)
  SBS35_detect(tt)
  
  cpu_fig_and_table()
  
  sup_table_s6(tt, "SBS")
  sup_table_s6(tt, "indel")
}  


# Calling plotting functions --------------------------------------------------
load(outpath("supplementary_table_s4.Rdata")) 
# Load supplementary_table_s4 computed in summarize_level1_dirs.R
all_figs_and_tables_this_file(supplementary_table_s4)
# 
