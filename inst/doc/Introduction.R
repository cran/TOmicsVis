## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
# 1. TOmicsVis
# install.packages("devtools")
# devtools::install_github("benben-miao/TOmicsVis")
library(TOmicsVis)

# 2. Extra package
# install.packages("ggplot2")
library(ggplot2)

## -----------------------------------------------------------------------------
# 1. Load deg_data example datasets
data(deg_data)

# 2. Run volcano_plot plot function
volcano_plot(
  deg_data,
  log2fc_cutoff = 1,
  pq_value = "pvalue",
  pq_cutoff = 0.005,
  cutoff_line = "longdash",
  point_shape = "large_circle",
  point_size = 1,
  point_alpha = 0.5,
  color_normal = "#888888",
  color_log2fc = "#008000",
  color_pvalue = "#0088ee",
  color_Log2fc_p = "#ff0000",
  label_size = 3,
  boxed_labels = FALSE,
  draw_connectors = FALSE,
  legend_pos = "right"
)

