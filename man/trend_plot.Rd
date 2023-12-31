% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trend_plot.R
\name{trend_plot}
\alias{trend_plot}
\title{Trend plot for visualizing gene expression trend profile in multiple traits.}
\usage{
trend_plot(
  data,
  scale_method = "centerObs",
  miss_value = "exclude",
  line_alpha = 0.5,
  show_points = TRUE,
  show_boxplot = TRUE,
  num_column = 1,
  xlab = "Traits",
  ylab = "Genes Expression",
  sci_fill_color = "Sci_AAAS",
  sci_fill_alpha = 0.8,
  sci_color_alpha = 0.8,
  legend_pos = "right",
  legend_dir = "vertical",
  ggTheme = "theme_light"
)
}
\arguments{
\item{data}{Dataframe: Shared degs of all paired comparisons in all groups expression dataframe of RNA-Seq. (1st-col: Genes, 2nd-col~n-1-col: Groups, n-col: Pathways).}

\item{scale_method}{Character: data scale methods. Default: "globalminmax" (global min and max values), options: "std" (standard), "robust", "uniminmax" (unique min and max values), "globalminmax", "center", "centerObs" (center observes).}

\item{miss_value}{Character: deal method for missing values. Default: "exclude", options: "exclude", "mean", "median", "min10", "random".}

\item{line_alpha}{Numeric: lines color alpha. Default: 0.50, min: 0.00, max: 1.00.}

\item{show_points}{Logical: show points at trait node. Default: TRUE, options: TRUE, FALSE.}

\item{show_boxplot}{Logical: show boxplot at trait node. Default: TRUE, options: TRUE, FALSE.}

\item{num_column}{Logical: column number. Default: 2, min: 1, max: null.}

\item{xlab}{Character: x label. Default: "Traits".}

\item{ylab}{Character: y label. Default: "Genes Expression".}

\item{sci_fill_color}{Character: ggsci color pallet. Default: "Sci_AAAS", options: "Sci_AAAS", "Sci_NPG", "Sci_Simpsons", "Sci_JAMA", "Sci_GSEA", "Sci_Lancet", "Sci_Futurama", "Sci_JCO", "Sci_NEJM", "Sci_IGV", "Sci_UCSC", "Sci_D3", "Sci_Material".}

\item{sci_fill_alpha}{Numeric: ggsci fill color alpha. Default: 0.50, min: 0.00, max: 1.00.}

\item{sci_color_alpha}{Numeric: ggsci border color alpha. Default: 1.00, min: 0.00, max: 1.00.}

\item{legend_pos}{Character: legend position. Default: "right", options: "none", "left", "right", "bottom", "top".}

\item{legend_dir}{Character: legend direction. Default: "vertical", options: "horizontal", "vertical".}

\item{ggTheme}{Character: ggplot2 themes. Default: "theme_light", options: "theme_default", "theme_bw", "theme_gray", "theme_light", "theme_linedraw", "theme_dark", "theme_minimal", "theme_classic", "theme_void"}
}
\value{
Plot: box plot support two levels and multiple groups with P value.
}
\description{
Trend plot for visualizing gene expression trend profile in multiple traits.
}
\examples{
# 1. Library TOmicsVis package
library(TOmicsVis)

# 2. Use example dataset
data(gene_expression3)
head(gene_expression3)

# 3. Default parameters
trend_plot(gene_expression3[1:50,])

# 4. Set line_alpha = 0.30
trend_plot(gene_expression3[1:50,], line_alpha = 0.30)

# 5. Set sci_fill_color = "Sci_NPG"
trend_plot(gene_expression3[1:50,], sci_fill_color = "Sci_NPG")

}
\author{
benben-miao
}
