## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  # Start shiny application.
#  TOmicsVis::tomicsvis()

## ----eval=FALSE---------------------------------------------------------------
#  # Install required packages from Bioconductor
#  install.packages("BiocManager")
#  BiocManager::install(c("ComplexHeatmap", "EnhancedVolcano", "clusterProfiler", "enrichplot", "impute", "preprocessCore", "Mfuzz"))

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("benben-miao/TOmicsVis")
#  
#  # Resolve network by GitClone
#  devtools::install_git("https://gitclone.com/github.com/benben-miao/TOmicsVis.git")

## ----eval=FALSE---------------------------------------------------------------
#  # Install from CRAN
#  install.packages("TOmicsVis")

## -----------------------------------------------------------------------------
# 1. Library TOmicsVis package
library(TOmicsVis)

# 2. Extra package
# install.packages("ggplot2")
library(ggplot2)

## ----quantile_plot, warning=FALSE, fig.width=10.00, fig.height=6.18-----------
# 1. Load example datasets
data(weight_sex)
head(weight_sex)

# 2. Run quantile_plot plot function
quantile_plot(
  data = weight_sex,
  my_shape = "fill_circle",
  point_size = 1.5,
  conf_int = TRUE,
  conf_level = 0.95,
  split_panel = "Split_Panel",
  legend_pos = "right",
  legend_dir = "vertical",
  sci_fill_color = "Sci_NPG",
  sci_color_alpha = 0.75,
  ggTheme = "theme_light"
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::quantile_plot

## ----box_plot, warning=FALSE, fig.width=10.00, fig.height=6.18----------------
# 1. Load example datasets
data(traits_sex)
head(traits_sex)

# 2. Run box_plot plot function
box_plot(
  data = traits_sex,
  test_method = "t.test",
  test_label = "p.format",
  notch = TRUE,
  group_level = "Three_Column",
  add_element = "jitter",
  my_shape = "fill_circle",
  sci_fill_color = "Sci_AAAS",
  sci_fill_alpha = 0.5,
  sci_color_alpha = 1,
  legend_pos = "right",
  legend_dir = "vertical",
  ggTheme = "theme_light"
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::box_plot

## ----violin_plot, warning=FALSE, fig.width=10.00, fig.height=6.18-------------
# 1. Load example datasets
data(traits_sex)

# 2. Run violin_plot plot function
violin_plot(
  data = traits_sex,
  test_method = "t.test",
  test_label = "p.format",
  group_level = "Three_Column",
  violin_orientation = "vertical",
  add_element = "boxplot",
  element_alpha = 0.5,
  my_shape = "plus_times",
  sci_fill_color = "Sci_AAAS",
  sci_fill_alpha = 0.5,
  sci_color_alpha = 1,
  legend_pos = "right",
  legend_dir = "vertical",
  ggTheme = "theme_light"
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::violin_plot

## ----survival_plot, warning=FALSE, fig.width=10.00, fig.height=6.18-----------
# 1. Load example datasets
data(survival_data)
head(survival_data)

# 2. Run survival_plot plot function
survival_plot(
  data = survival_data,
  curve_function = "pct",
  conf_inter = TRUE,
  interval_style = "ribbon",
  risk_table = TRUE,
  num_censor = TRUE,
  sci_palette = "aaas",
  ggTheme = "theme_light",
  x_start = 0,
  y_start = 0,
  y_end = 100,
  x_break = 10,
  y_break = 10
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::survival_plot

## ----corr_heatmap, fig.width=10.00, fig.height=6.18---------------------------
# 1. Load example dataset
data(gene_expression)
head(gene_expression)

# 2. Run corr_heatmap plot function
corr_heatmap(
  data = gene_expression,
  corr_method = "pearson",
  cell_shape = "square",
  fill_type = "full",
  lable_size = 3,
  axis_angle = 45,
  axis_size = 12,
  lable_digits = 3,
  color_low = "blue",
  color_mid = "white",
  color_high = "red",
  outline_color = "white",
  ggTheme = "theme_light"
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::corr_heatmap

## ----pca_analysis, fig.width=10.00, fig.height=6.18---------------------------
# 1. Load example datasets
data(gene_expression)

data(samples_groups)
head(samples_groups)

# 2. Run pca_analysis plot function
res <- pca_analysis(gene_expression, samples_groups)
head(res)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::pca_analysis

## ----pca_plot, fig.width=10.00, fig.height=6.18-------------------------------
# 1. Load example datasets
data(gene_expression)

data(samples_groups)
head(samples_groups)

# 2. Run pca_plot plot function
pca_plot(
  sample_gene = gene_expression,
  group_sample = samples_groups,
  xPC = 1,
  yPC = 2,
  point_size = 5,
  text_size = 5,
  fill_alpha = 0.10,
  border_alpha = 0.00,
  legend_pos = "right",
  legend_dir = "vertical",
  ggTheme = "theme_light"
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::pca_plot

## ----tsne_analysis, warning=FALSE, fig.width=10.00, fig.height=6.18-----------
# 1. Load example datasets
data(gene_expression)
data(samples_groups)

# 2. Run tsne_analysis plot function
res <- tsne_analysis(gene_expression, samples_groups)
head(res)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::tsne_analysis

## ----tsne_plot, warning=FALSE, fig.width=10.00, fig.height=6.18---------------
# 1. Load example datasets
data(gene_expression)
data(samples_groups)

# 2. Run tsne_plot plot function
tsne_plot(
  sample_gene = gene_expression,
  group_sample = samples_groups,
  seed = 1,
  multi_shape = FALSE,
  point_size = 5,
  point_alpha = 0.8,
  text_size = 5,
  text_alpha = 0.80,
  fill_alpha = 0.10,
  border_alpha = 0.00,
  sci_fill_color = "Sci_AAAS",
  legend_pos = "right",
  legend_dir = "vertical",
  ggTheme = "theme_light"
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::tsne_plot

## ----umap_analysis, warning=FALSE, fig.width=10.00, fig.height=6.18-----------
# 1. Load example datasets
data(gene_expression)
data(samples_groups)

# 2. Run tsne_plot plot function
res <- umap_analysis(gene_expression, samples_groups)
head(res)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::umap_analysis

## ----umap_plot, warning=FALSE, fig.width=10.00, fig.height=6.18---------------
# 1. Load example datasets
data(gene_expression)
data(samples_groups)

# 2. Run tsne_plot plot function
umap_plot(
  sample_gene = gene_expression,
  group_sample = samples_groups,
  seed = 1,
  multi_shape = TRUE,
  point_size = 5,
  point_alpha = 1,
  text_size = 5,
  text_alpha = 0.80,
  fill_alpha = 0.00,
  border_alpha = 0.00,
  sci_fill_color = "Sci_AAAS",
  legend_pos = "right",
  legend_dir = "vertical",
  ggTheme = "theme_light"
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::umap_plot

## ----dendro_plot, fig.width=10.00, fig.height=6.18----------------------------
# 1. Load example datasets
data(gene_expression)

# 2. Run plot function
dendro_plot(
  data = gene_expression,
  dist_method = "euclidean",
  hc_method = "ward.D2",
  tree_type = "rectangle",
  k_num = 5,
  palette = "npg",
  color_labels_by_k = TRUE,
  horiz = FALSE,
  label_size = 1,
  line_width = 1,
  rect = TRUE,
  rect_fill = TRUE,
  xlab = "Samples",
  ylab = "Height",
  ggTheme = "theme_light"
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::dendro_plot

## ----venn_plot, warning=FALSE, fig.width=10.00, fig.height=6.18---------------
# 1. Load example datasets
data(degs_lists)
head(degs_lists)

# 2. Run venn_plot plot function
venn_plot(
  data = degs_lists,
	title_size = 1,
	label_show = TRUE,
	label_size = 0.8,
	border_show = TRUE,
	line_type = "longdash",
	ellipse_shape = "circle",
	sci_fill_color = "Sci_AAAS",
	sci_fill_alpha = 0.65
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::venn_plot

## ----upsetr_plot, warning=FALSE, fig.width=10.00, fig.height=6.18-------------
# 1. Load example datasets
data(degs_lists)
head(degs_lists)

# 2. Run upsetr_plot plot function
upsetr_plot(
  data = degs_lists,
  sets_num = 4,
  keep_order = FALSE,
  order_by = "freq",
  decrease = TRUE,
  mainbar_color = "#006600",
  number_angle = 45,
  matrix_color = "#cc0000",
  point_size = 4.5,
  point_alpha = 0.5,
  line_size = 0.8,
  shade_color = "#cdcdcd",
  shade_alpha = 0.5,
  setsbar_color = "#000066",
  setsnum_size = 6,
  text_scale = 1.2
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::upsetr_plot

## ----flower_plot, warning=FALSE, fig.width=10.00, fig.height=6.18-------------
# 1. Load example datasets
data(degs_lists)

# 2. Run plot function
flower_plot(
  flower_dat = degs_lists,
  angle = 90,
  a = 1,
  b = 2,
  r = 1,
  ellipse_col_pal = "Spectral",
  circle_col = "white",
  label_text_cex = 1
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::flower_plot

## ----volcano_plot, warning=FALSE, fig.width=10.00, fig.height=6.18------------
# 1. Load example datasets
data(degs_stats)
head(degs_stats)

# 2. Run volcano_plot plot function
volcano_plot(
  data = degs_stats,
  title = "CT-vs-LT12",
  log2fc_cutoff = 1,
  pq_value = "pvalue",
  pq_cutoff = 0.05,
  cutoff_line = "longdash",
  point_shape = "large_circle",
  point_size = 2,
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

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::volcano_plot

## ----ma_plot, warning=FALSE, fig.width=10.00, fig.height=6.18-----------------
# 1. Load example datasets
data(degs_stats2)
head(degs_stats2)

# 2. Run volcano_plot plot function
ma_plot(
  data = degs_stats2,
  foldchange = 2,
  fdr_value = 0.05,
  point_size = 3.0,
  color_up = "#FF0000",
  color_down = "#008800",
  color_alpha = 0.5,
  top_method = "fc",
  top_num = 20,
  label_size = 8,
  label_box = TRUE,
  title = "CT-vs-LT12",
  xlab = "Log2 mean expression",
  ylab = "Log2 fold change",
  ggTheme = "theme_light"
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::ma_plot

## ----heatmap_group, warning=FALSE, fig.width=10.00, fig.height=6.18-----------
# 1. Load example datasets
data(gene_expression2)
data(samples_groups)

# 2. Run heatmap_group plot function
heatmap_group(
  sample_gene = gene_expression2[1:30,],
  group_sample = samples_groups,
  scale_data = "row",
  clust_method = "complete",
  border_show = TRUE,
  border_color = "#ffffff",
  value_show = TRUE,
  value_decimal = 2,
  value_size = 5,
  axis_size = 8,
  cell_height = 10,
  low_color = "#00880055",
  mid_color = "#ffffff",
  high_color = "#ff000055",
  na_color = "#ff8800",
  x_angle = 45
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::heatmap_group

## ----circos_heatmap, warning=FALSE, fig.width=10.00, fig.height=6.18----------
# 1. Load example datasets
data(gene_expression2)
head(gene_expression2)

# 2. Run circos_heatmap plot function
circos_heatmap(
  data = gene_expression2[1:50,],
  low_color = "#0000ff",
  mid_color = "#ffffff",
  high_color = "#ff0000",
  gap_size = 25,
  cluster_run = TRUE,
  cluster_method = "complete",
  distance_method = "euclidean",
  dend_show = "inside",
  dend_height = 0.2,
  track_height = 0.3,
  rowname_show = "outside",
  rowname_size = 0.8
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::circos_heatmap

## ----chord_plot, warning=FALSE, fig.width=10.00, fig.height=6.18--------------
# 1. Load chord_data example datasets
data(gene_expression2)
head(gene_expression2)

# 2. Run chord_plot plot function
chord_plot(
  data = gene_expression2[1:30,],
  multi_colors = "VividColors",
  color_seed = 10,
  color_alpha = 0.3,
  link_visible = TRUE,
  link_dir = -1,
  link_type = "diffHeight",
  sector_scale = "Origin",
  width_circle = 3,
  dist_name = 3,
  label_dir = "Vertical",
  dist_label = 0.3,
  label_scale = 0.8
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::chord_plot

## ----gene_rank_plot, warning=FALSE, fig.width=10.00, fig.height=6.18----------
# 1. Load example datasets
data(degs_stats)

# 2. Run plot function
gene_rank_plot(
  data = degs_stats,
  log2fc = 1,
  palette = "Spectral",
  top_n = 10,
  genes_to_label = NULL,
  label_size = 5,
  base_size = 12,
  title = "Gene ranking dotplot",
  xlab = "Ranking of differentially expressed genes",
  ylab = "Log2FoldChange"
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::gene_rank_plot

## ----gene_cluster_trend, warning=FALSE, fig.width=10.00, fig.height=6.18------
# 1. Load example datasets
data(gene_expression3)

# 2. Run plot function
gene_cluster_trend(
  data = gene_expression3[,-7],
  thres = 0.25,
  min_std = 0.2,
  palette = "PiYG",
  cluster_num = 4
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::gene_cluster_trend

## ----trend_plot, warning=FALSE, fig.width=10.00, fig.height=6.18--------------
# 1. Load example datasets
data(gene_expression3)
head(gene_expression3)

# 2. Run trend_plot plot function
trend_plot(
  data = gene_expression3[1:100,],
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

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::trend_plot

## ----wgcna_pipeline, warning=FALSE, fig.width=10.00, fig.height=6.18----------
# 1. Load wgcna_pipeline example datasets
data(gene_expression)
head(gene_expression)

data(samples_groups)
head(samples_groups)

# 2. Run wgcna_pipeline plot function
# wgcna_pipeline(gene_expression[1:3000,], samples_groups)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::wgcna_pipeline

## ----network_plot, warning=FALSE, fig.width=10.00, fig.height=6.18------------
# 1. Load example datasets
data(network_data)
head(network_data)

# 2. Run network_plot plot function
network_plot(
  data = network_data,
  calc_by = "degree",
  degree_value = 0.5,
  normal_color = "#008888cc",
  border_color = "#FFFFFF",
  from_color = "#FF0000cc",
  to_color = "#008800cc",
  normal_shape = "circle",
  spatial_shape = "circle",
  node_size = 25,
  lable_color = "#FFFFFF",
  label_size = 0.5,
  edge_color = "#888888",
  edge_width = 1.5,
  edge_curved = TRUE,
  net_layout = "layout_on_sphere"
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::network_plot

## ----heatmap_cluster, warning=FALSE, fig.width=10.00, fig.height=6.18---------
# 1. Load example datasets
data(gene_expression2)
head(gene_expression2)

# 2. Run network_plot plot function
heatmap_cluster(
  data = gene_expression2,
  dist_method = "euclidean",
  hc_method = "average",
  k_num = 5,
  show_rownames = FALSE,
  palette = "RdBu",
  cluster_pal = "Set1",
  border_color = "#ffffff",
  angle_col = 45,
  label_size = 10,
  base_size = 12,
  line_color = "#0000cd",
  line_alpha = 0.2,
  summary_color = "#0000cd",
  summary_alpha = 0.8
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::heatmap_cluster

## ----go_enrich, warning=FALSE, fig.width=10.00, fig.height=6.18---------------
# 1. Load example datasets
data(gene_go_kegg)
head(gene_go_kegg)

# 2. Run go_enrich analysis function
res <- go_enrich(
  go_anno = gene_go_kegg[,-5],
  degs_list = gene_go_kegg[100:200,1],
  padjust_method = "fdr",
  pvalue_cutoff = 0.05,
  qvalue_cutoff = 0.05
)
head(res)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::go_enrich

## ----go_enrich_stat, warning=FALSE, fig.width=10.00, fig.height=6.18----------
# 1. Load example datasets
data(gene_go_kegg)

# 2. Run go_enrich_stat analysis function
go_enrich_stat(
  go_anno = gene_go_kegg[,-5],
  degs_list = gene_go_kegg[100:200,1],
  padjust_method = "fdr",
  pvalue_cutoff = 0.05,
  qvalue_cutoff = 0.05,
  max_go_item = 15,
  strip_fill = "#CDCDCD",
  xtext_angle = 45,
  sci_fill_color = "Sci_AAAS",
  sci_fill_alpha = 0.8,
  ggTheme = "theme_light"
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::go_enrich_stat

## ----go_enrich_bar, warning=FALSE, fig.width=10.00, fig.height=6.18-----------
# 1. Load example datasets
data(gene_go_kegg)

# 2. Run go_enrich_bar analysis function
go_enrich_bar(
  go_anno = gene_go_kegg[,-5],
  degs_list = gene_go_kegg[100:200,1],
  padjust_method = "fdr",
  pvalue_cutoff = 0.05,
  qvalue_cutoff = 0.05,
  sign_by = "p.adjust",
  category_num = 30,
  font_size = 12,
  low_color = "#ff0000aa",
  high_color = "#008800aa",
  ggTheme = "theme_light"
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::go_enrich_bar

## ----go_enrich_dot, warning=FALSE, fig.width=10.00, fig.height=6.18-----------
# 1. Load example datasets
data(gene_go_kegg)

# 2. Run go_enrich_dot analysis function
go_enrich_dot(
  go_anno = gene_go_kegg[,-5],
  degs_list = gene_go_kegg[100:200,1],
  padjust_method = "fdr",
  pvalue_cutoff = 0.05,
  qvalue_cutoff = 0.05,
  sign_by = "p.adjust",
  category_num = 30,
  font_size = 12,
  low_color = "#ff0000aa",
  high_color = "#008800aa",
  ggTheme = "theme_light"
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::go_enrich_dot

## ----go_enrich_net, warning=FALSE, fig.width=10.00, fig.height=6.18-----------
# 1. Load example datasets
data(gene_go_kegg)

# 2. Run go_enrich_net analysis function
go_enrich_net(
  go_anno = gene_go_kegg[,-5],
  degs_list = gene_go_kegg[100:200,1],
  padjust_method = "fdr",
  pvalue_cutoff = 0.05,
  qvalue_cutoff = 0.05,
  category_num = 20,
  net_layout = "circle",
  net_circular = TRUE,
  low_color = "#ff0000aa",
  high_color = "#008800aa"
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::go_enrich_net

## ----kegg_enrich, warning=FALSE, fig.width=10.00, fig.height=6.18-------------
# 1. Load example datasets
data(gene_go_kegg)
head(gene_go_kegg)

# 2. Run go_enrich analysis function
res <- kegg_enrich(
  kegg_anno = gene_go_kegg[,c(1,5)],
  degs_list = gene_go_kegg[100:200,1],
  padjust_method = "fdr",
  pvalue_cutoff = 0.05,
  qvalue_cutoff = 0.05
)
head(res)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::kegg_enrich

## ----kegg_enrich_bar, warning=FALSE, fig.width=10.00, fig.height=6.18---------
# 1. Load example datasets
data(gene_go_kegg)

# 2. Run kegg_enrich_bar analysis function
kegg_enrich_bar(
  kegg_anno = gene_go_kegg[,c(1,5)],
  degs_list = gene_go_kegg[100:200,1],
  padjust_method = "fdr",
  pvalue_cutoff = 0.05,
  qvalue_cutoff = 0.05,
  sign_by = "p.adjust",
  category_num = 30,
  font_size = 12,
  low_color = "#ff0000aa",
  high_color = "#008800aa",
  ggTheme = "theme_light"
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::kegg_enrich_bar

## ----kegg_enrich_dot, warning=FALSE, fig.width=10.00, fig.height=6.18---------
# 1. Load example datasets
data(gene_go_kegg)

# 2. Run kegg_enrich_dot analysis function
kegg_enrich_dot(
  kegg_anno = gene_go_kegg[,c(1,5)],
  degs_list = gene_go_kegg[100:200,1],
  padjust_method = "fdr",
  pvalue_cutoff = 0.05,
  qvalue_cutoff = 0.05,
  sign_by = "p.adjust",
  category_num = 30,
  font_size = 12,
  low_color = "#ff0000aa",
  high_color = "#008800aa",
  ggTheme = "theme_light"
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::kegg_enrich_dot

## ----kegg_enrich_net, warning=FALSE, fig.width=10.00, fig.height=6.18---------
# 1. Load example datasets
data(gene_go_kegg)

# 2. Run kegg_enrich_net analysis function
kegg_enrich_net(
  kegg_anno = gene_go_kegg[,c(1,5)],
  degs_list = gene_go_kegg[100:200,1],
  padjust_method = "fdr",
  pvalue_cutoff = 0.05,
  qvalue_cutoff = 0.05,
  category_num = 20,
  net_layout = "circle",
  net_circular = TRUE,
  low_color = "#ff0000aa",
  high_color = "#008800aa"
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::kegg_enrich_net

## ----table_split, warning=FALSE, fig.width=10.00, fig.height=6.18-------------
# 1. Load example datasets
data(gene_go_kegg2)
head(gene_go_kegg2)

# 2. Run table_split function
res <- table_split(
  data = gene_go_kegg2,
  grouped_var = "go_category",
  value_var = "go_term",
  miss_drop = TRUE
)
head(res)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::table_split

## ----table_merge, warning=FALSE, fig.width=10.00, fig.height=6.18-------------
# 1. Load example datasets
data(gene_go_kegg)
head(gene_go_kegg)

# 2. Run function
res <- table_merge(
  data = gene_go_kegg,
  merge_vars = c("biological_process", "cellular_component", "molecular_function"),
  new_var = "go_category",
  new_value = "go_term",
  na_remove = FALSE
)
head(res)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::table_merge

## ----table_filter, warning=FALSE, fig.width=10.00, fig.height=6.18------------
# 1. Load example datasets
data(traits_sex)
head(traits_sex)

# 2. Run function
res <- table_filter(
	data = traits_sex, 
	Sex == "Male" & Traits == "Weight" & Value > 40
	)
head(res)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::table_filter

## ----table_cross, warning=FALSE, fig.width=10.00, fig.height=6.18-------------
# 1. Load example datasets
data(gene_expression2)
head(gene_expression2)

data(gene_go_kegg)
head(gene_go_kegg)

# 2. Run function
res <- table_cross(
  data1 = gene_expression2,
  data2 = gene_go_kegg,
  inter_var = "Genes",
  left_index = TRUE,
  right_index = TRUE
)
head(res)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::table_cross

