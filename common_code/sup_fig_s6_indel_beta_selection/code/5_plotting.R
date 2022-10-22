# Please run this script from the top directory
if (basename(getwd()) != "Liu_et_al_Sup_Files") {
  stop("Please run from top level directory, Liu_et_al_Sup_Files")
}

# 0. Install and load dependencies --------------------------------------------

require(ggplot2)
require(dplyr)
require(ggbeeswarm)
require(ggpubr)
require(data.table)
require(viridis)

# 1. Specify global variables -------------------------------------------------

home_for_summary <- "other_analyses/indel_beta_selection/summary/"



# 2. Define plotting function -------------------------------------------------

# A plotting function which accepts a data frame sliced
# form all_results.csv as input,
# and a ggplot object contains beeswarm plot for a variable as output.
#
# DF - Input data.frame object.
#
# var.name - Name of variable (as a character) to be plotted (y).
# This variable will be plotted against noise level (x)
# and computational approaches (color and shape)
#
# var.title - The title to display on y axis.
#
# inputRange - Range of y axis for variable. Passed to limits argument of
# scale_y_continuous().
#
# Sci - Whether to enable scientific notation.
plotFunc <- function(DF, var.name, var.title, inputRange, sci = FALSE){
  
  DF1 <- DF[!is.na(DF[[var.name]]),]
  
  ggObj <- ggplot2::ggplot(
    data = DF1,
    mapping = aes_string(x = "beta_level",
                         y = var.name)) +
    # Add a dot plot
    # Dot plot is a better visualization
    # compared to boxplot or violin plot,
    # as there are only 20 points in each group
    ggbeeswarm::geom_beeswarm(
      # All dots have shape "open circle",
      # to be consistent with top-level summary plot.
      shape = 19,
      # Shrink the size of dots.
      size = 0.7,
      # Make dots in each aesthetic group
      # dodge from one another.
      # dodge.width equals to the distance between
      # the leftmost and the rightmost point
      # in each group / the distance between two group
      # tick marks.
      dodge.width = 0.6,
      # The jitter is on the X axis rather than on
      # Y axis.
      # This can keep the value of composite measure.
      groupOnX = T) +
    # Order the groups by variable "beta_level".
    ggplot2::scale_x_discrete(limits = DF$beta_level %>% 
                                unique() %>% 
                                gtools::mixedsort()) +
    # Change the range for different scatterplot to be the same.
    ggplot2::scale_y_continuous(
      labels = function(x) format(x, scientific = sci),
      limits = inputRange) +
    # Scale color by viridis
    viridis::scale_color_viridis(
      # Forward direction
      direction = 1,
      # Alpha transparency
      # 1: atransparent
      alpha = 1,
      # Preventing the two extremities
      begin = 0.1,
      end = 0.9,
      discrete = T,
      # Palette version
      option = "H") +
    # Change axis and legend titles
    ggplot2::labs(x = "beta",
                  y = var.title,
                  color = "Computational Approach") +
    # Rotate axis.text.x 90 degrees,
    # move axis.text.x right below the tick marks,
    # and remove legends.
    ggplot2::theme(
      # Change fill and border color of background rectangle
      panel.background = ggplot2::element_rect(
        # Fill color
        fill = "white", 
        # Border color
        colour = "black"),
      # Remove white grid lines
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        # Rotate the axis.text.x
        angle = 30,
        # move axis.text.x right below the tick marks
        hjust = 1,vjust = 1),
      # Make font size of facet label smaller.
      strip.text = ggplot2::element_text(size = 10),
      # remove legends.
      legend.position = "top") +
    # Let legend span across two rows.
    ggplot2::guides(
      color = guide_legend(nrow = 2))
  
  return(ggObj)
}



# 3. Plot extraction measures -------------------------------------------------

# a. Import results of extraction measures.
file <- paste0(home_for_summary, "/all_results.csv")
DF <- read.csv(file, header = T)


# b. Data pre-processing.
#
# Plot for PPV, TPR, mean Cosine similarity,
# and Composite Measure equals to the sum of
# these 3 measures.
# Calculate composite measure
DF <- DF %>% mutate(Composite = PPV + TPR + aver_Sim_TP_only)

# Change tool names to values of beta argument
DF <- DF %>% mutate(beta_level =  gsub("mSigHdp.beta", "", Approach),
                    .keep = "unused",
                    .before = "Noise_level")


# c. Plotting of 4 panels, and saving arranged plot.
var.names <- c("aver_Sim_TP_only","PPV","TPR")
axis.titles <- c("Composite" = "Composite measure",
                 "aver_Sim_TP_only" = "Mean cosine similarity",
                 "PPV" = "PPV",
                 "TPR" = "TPR")

figs <- list()
for(vn in c("Composite", var.names)){
  # Round range of to 1 decimal place.
  range.vals <- DF %>% select(all_of(vn)) %>%
    unlist() %>% range()
    if (vn != "aver_Sim_TP_only") {
      range.vals <-
        c(floor(10 * range.vals[1])/10, ceiling(10 * range.vals[2])/10)
    } else {
      range.vals <-
        c(floor(100 * range.vals[1])/100, ceiling(100 * range.vals[2])/100)
    }
  # Generate ggplot object
  figs[[vn]] <- plotFunc(DF, vn, axis.titles[vn], range.vals)
}
# Arrange and save the plot
arr.figs <- ggpubr::ggarrange(plotlist = figs, common.legend = T)
ggplot2::ggsave(
  filename = paste0(home_for_summary,"/ID.testing.diff.beta.with.SP.pdf"),
  plot = arr.figs,
  width = unit(9, "inch"),
  height = unit(6, "inch"))


# d. Plotting of 4 panels, and saving arranged plot,
# yet with SigProfilerExtractor removed. 
DF_no_SP <- DF %>% filter(beta_level != "SigProfilerExtractor")

figs <- list()
for(vn in c("Composite", var.names)){
  # Round range of to 1 decimal place.
  range.vals <- DF_no_SP %>% select(all_of(vn)) %>%
    unlist() %>% range()
  if (vn != "aver_Sim_TP_only") {
    range.vals <-
      c(floor(10 * range.vals[1])/10, ceiling(10 * range.vals[2])/10)
  } else {
    range.vals <-
      c(floor(100 * range.vals[1])/100, ceiling(100 * range.vals[2])/100)
  }
  # Generate ggplot object
  figs[[vn]] <- plotFunc(DF_no_SP, vn, axis.titles[vn], range.vals)
}
# Arrange and save the plot
arr.figs <- ggpubr::ggarrange(plotlist = figs, common.legend = T)
arr.figs <- ggpubr::annotate_figure(arr.figs, fig.lab = "burnin = 30k", fig.lab.face = "bold")
ggplot2::ggsave(
  filename = paste0(home_for_summary,"/ID.testing.diff.beta.pdf"),
  plot = arr.figs,
  width = unit(9, "inch"),
  height = unit(6, "inch"))
