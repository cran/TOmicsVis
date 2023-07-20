## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

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
# 1. Load box_data example datasets
data(quantile_data)

# 2. Run quantile_plot plot function
quantile_plot(
  quantile_data,
  my_shape = "fill_circle",
  point_size = 1.5,
  conf_int = TRUE,
  conf_level = 0.95,
  split_panel = "One_Panel",
  legend_pos = "right",
  legend_dir = "vertical",
  sci_fill_color = "Sci_AAAS",
  sci_color_alpha = 0.75,
  ggTheme = "theme_light"
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::quantile_plot

## ----corr_heatmap, fig.width=10.00, fig.height=6.18---------------------------
# 1. Load gene_exp example dataset
data(gene_exp)
head(gene_exp)

# 2. Run corr_heatmap plot function
corr_heatmap(
  gene_exp,
  corr_method = "pearson",
  cell_shape = "square",
  fill_type = "full",
  lable_size = 3,
  lable_digits = 3,
  color_low = "blue",
  color_mid = "white",
  color_high = "red",
  ggTheme = "theme_light"
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::corr_heatmap

## ----pca_plot, fig.width=10.00, fig.height=6.18-------------------------------
# 1. Load pca_sample_gene and pca_group_sample example datasets
data(pca_sample_gene)
data(pca_group_sample)

# 2. Run pca_plot plot function
pca_plot(
  pca_sample_gene,
  pca_group_sample,
  point_size = 5,
  text_size = 5,
  ellipse_alpha = 0.3,
  legend_pos = "right",
  legend_dir = "vertical",
  ggTheme = "theme_light"
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::pca_plot

## ----dendro_plot, fig.width=10.00, fig.height=6.18----------------------------
# 1. Load example datasets
data(gene_exp)

# 2. Run plot function
dendro_plot(
  gene_exp,
  dist_method = "euclidean",
  hc_method = "average",
  tree_type = "rectangle",
  k_num = 3,
  palette = "npg",
  color_labels_by_k = TRUE,
  horiz = TRUE,
  label_size = 0.8,
  line_width = 0.7,
  rect = TRUE,
  rect_fill = TRUE,
  title = "Cluster Dendrogram",
  xlab = "",
  ylab = "Height"
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::dendro_plot

## ----box_plot, warning=FALSE, fig.width=10.00, fig.height=6.18----------------
# 1. Load box_data example datasets
data(box_data)

# 2. Run box_plot plot function
box_plot(
  box_data,
  test_method = "t.test",
  test_label = "p.format",
  notch = TRUE,
  group_level = "Three_Column",
  add_element = "dotplot",
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
# 1. Load box_data example datasets
data(box_data)

# 2. Run violin_plot plot function
violin_plot(
  box_data,
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
# 1. Load survival_plot example datasets
data(survival_data)

# 2. Run survival_plot plot function
survival_plot(
  survival_data,
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
  x_break = 100,
  y_break = 25
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::survival_plot

## ----tsne_plot, warning=FALSE, fig.width=10.00, fig.height=6.18---------------
# 1. Load tsne_plot example datasets
data(tsne_data)

# 2. Run tsne_plot plot function
tsne_plot(
  tsne_data,
  seed = 5,
  point_size = 4,
  point_alpha = 0.8,
  text_size = 2,
  text_alpha = 0.8,
  ci_level = 0.95,
  ellipse_alpha = 0.3,
  sci_fill_color = "Sci_JAMA",
  sci_color_alpha = 0.9,
  legend_pos = "right",
  legend_dir = "vertical",
  ggTheme = "theme_light"
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::tsne_plot

## ----venn_plot, warning=FALSE, fig.width=10.00, fig.height=6.18---------------
# 1. Load venn_data example datasets
data(venn_data)

# 2. Run venn_plot plot function
venn_plot(
  venn_data,
  line_type = "blank",
  ellipse_shape = "circle",
  sci_fill_color = "Sci_AAAS",
  sci_fill_alpha = 0.65
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::venn_plot

## ----flower_plot, warning=FALSE, fig.width=10.00, fig.height=6.18-------------
# 1. Load example datasets
data(venn_data)

# 2. Run plot function
flower_plot(
  venn_data,
  angle = 90,
  a = 0.5,
  b = 2,
  r = 1,
  ellipse_col_pal = "Spectral",
  circle_col = "white",
  label_text_cex = 1
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::flower_plot

## ----circos_heatmap, warning=FALSE, fig.width=10.00, fig.height=6.18----------
# 1. Load circos_heatmap_data example datasets
data(circos_heatmap_data)

# 2. Run circos_heatmap plot function
circos_heatmap(
  circos_heatmap_data,
  low_color = "#0000ff",
  mid_color = "#ffffff",
  high_color = "#ff0000",
  gap_size = 10,
  cluster_method = "complete",
  distance_method = "euclidean",
  dend_height = 0.2,
  rowname_size = 0.8
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::circos_heatmap

## ----volcano_plot, warning=FALSE, fig.width=10.00, fig.height=6.18------------
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

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::volcano_plot

## ----ma_plot, warning=FALSE, fig.width=10.00, fig.height=6.18-----------------
# 1. Load deg_data example datasets
data(deg_data2)

# 2. Run volcano_plot plot function
ma_plot(
  deg_data2,
  foldchange = 2,
  fdr_value = 0.05,
  point_size = 0.5,
  color_up = "#FF0000",
  color_down = "#008800",
  color_alpha = 0.5,
  top_method = "fc",
  top_num = 20,
  label_size = 8,
  label_box = TRUE,
  title = "Group1 -versus- Group2",
  xlab = "Log2 mean expression",
  ylab = "Log2 fold change",
  ggTheme = "theme_minimal"
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::ma_plot

## ----heatmap_group, warning=FALSE, fig.width=10.00, fig.height=6.18-----------
# 1. Load example datasets
data(heatmap_group_data)
head(heatmap_group_data)

# 2. Run heatmap_group plot function
heatmap_group(
  data = heatmap_group_data,
  scale_data = "none",
  clust_method = "complete",
  border_show = TRUE,
  value_show = TRUE,
  low_color = "#00880088",
  mid_color = "#ffffff",
  high_color = "#ff000088",
  na_color = "#ff8800",
  x_angle = 45
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::heatmap_group

## ----trend_plot, warning=FALSE, fig.width=10.00, fig.height=6.18--------------
# 1. Load chord_data example datasets
data(trend_data)

# 2. Run trend_plot plot function
trend_plot(
  trend_data,
  scale_method = "globalminmax",
  miss_value = "exclude",
  line_alpha = 0.5,
  show_points = TRUE,
  show_boxplot = TRUE,
  num_column = 2,
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

## ----gene_cluster_trend, warning=FALSE, fig.width=10.00, fig.height=6.18------
# 1. Load example datasets
data(gene_cluster_data)

# 2. Run plot function
gene_cluster_trend(
  gene_cluster_data,
  thres = 0.25,
  min_std = 0.2,
  palette = "PiYG",
  cluster_num = 4
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::gene_cluster_trend

## ----gene_rank_plot, warning=FALSE, fig.width=10.00, fig.height=6.18----------
# 1. Load example datasets
data(deg_data)

# 2. Run plot function
gene_rank_plot(
  data = deg_data,
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

## ----wgcna_pipeline, warning=FALSE, fig.width=10.00, fig.height=6.18----------
# 1. Load wgcna_pipeline example datasets
data(wgcna_gene_exp)
data(wgcna_sample_group)

# 2. Run wgcna_pipeline plot function
# wgcna_pipeline(wgcna_gene_exp, wgcna_sample_group)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::wgcna_pipeline

## ----network_plot, warning=FALSE, fig.width=10.00, fig.height=6.18------------
# 1. Load example datasets
data(network_data)
head(network_data)

# 2. Run network_plot plot function
network_plot(
  network_data,
  calcBy = "degree",
  degreeValue = 0.05,
  nodeColorNormal = "#00888888",
  nodeBorderColor = "#FFFFFF",
  nodeColorFrom = "#FF000088",
  nodeColorTo = "#00880088",
  nodeShapeNormal = "circle",
  nodeShapeSpatial = "csquare",
  nodeSize = 10,
  labelSize = 0.5,
  edgeCurved = TRUE,
  netLayout = "layout_on_sphere"
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::network_plot

## ----heatmap_cluster, warning=FALSE, fig.width=10.00, fig.height=6.18---------
# 1. Load example datasets
data(gene_exp2)
head(gene_exp2)

# 2. Run network_plot plot function
heatmap_cluster(
  data = gene_exp2,
  dist_method = "euclidean",
  hc_method = "average",
  k_num = 5,
  show_rownames = FALSE,
  palette = "Spectral",
  cluster_pal = "Set1",
  gaps_col = NULL,
  angle_col = 45,
  label_size = 10,
  base_size = 12
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::heatmap_cluster

## ----chord_plot, warning=FALSE, fig.width=10.00, fig.height=6.18--------------
# 1. Load chord_data example datasets
data(chord_data)

# 2. Run chord_plot plot function
chord_plot(
  chord_data,
  multi_colors = "RainbowColors",
  color_alpha = 0.5,
  link_visible = TRUE,
  link_dir = -1,
  link_type = "diffHeight",
  sector_scale = "Origin",
  width_circle = 3,
  dist_name = 3,
  label_dir = "Vertical",
  dist_label = 0.3
)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::chord_plot

## ----go_enrich, warning=FALSE, fig.width=10.00, fig.height=6.18---------------
# 1. Load example datasets
data(go_anno)
head(go_anno)

data(go_deg_fc)
head(go_deg_fc)

# 2. Run go_enrich analysis function
res <- go_enrich(
  go_anno,
  go_deg_fc,
  padjust_method = "fdr",
  pvalue_cutoff = 0.5,
  qvalue_cutoff = 0.5
)
head(res)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::go_enrich

## ----go_enrich_stat, warning=FALSE, fig.width=10.00, fig.height=6.18----------
# 1. Load example datasets
data(go_anno)
# head(go_anno)

data(go_deg_fc)
# head(go_deg_fc)

# 2. Run go_enrich_stat analysis function
go_enrich_stat(
  go_anno,
  go_deg_fc,
  padjust_method = "fdr",
  pvalue_cutoff = 0.5,
  qvalue_cutoff = 0.5,
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
data(go_anno)
# head(go_anno)

data(go_deg_fc)
# head(go_deg_fc)

# 2. Run go_enrich_bar analysis function
go_enrich_bar(
  go_anno,
  go_deg_fc,
  padjust_method = "fdr",
  pvalue_cutoff = 0.5,
  qvalue_cutoff = 0.5,
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
data(go_anno)
# head(go_anno)

data(go_deg_fc)
# head(go_deg_fc)

# 2. Run go_enrich_dot analysis function
go_enrich_dot(
  go_anno,
  go_deg_fc,
  padjust_method = "fdr",
  pvalue_cutoff = 0.5,
  qvalue_cutoff = 0.5,
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

## ----go_enrich_tree, warning=FALSE, fig.width=10.00, fig.height=6.18, eval=FALSE----
#  # 1. Load example datasets
#  data(go_anno)
#  # head(go_anno)
#  
#  data(go_deg_fc)
#  # head(go_deg_fc)
#  
#  # 2. Run go_enrich_tree analysis function
#  go_enrich_tree(
#    go_anno,
#    go_deg_fc,
#    padjust_method = "fdr",
#    pvalue_cutoff = 0.5,
#    qvalue_cutoff = 0.5,
#    sign_by = "p.adjust",
#    category_num = 20,
#    font_size = 4,
#    low_color = "#ff0000aa",
#    high_color = "#008800aa",
#    hclust_method = "complete",
#    ggTheme = "theme_light"
#  )

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::go_enrich_tree

## ----go_enrich_net, warning=FALSE, fig.width=10.00, fig.height=6.18-----------
# 1. Load example datasets
data(go_anno)
# head(go_anno)

data(go_deg_fc)
# head(go_deg_fc)

# 2. Run go_enrich_net analysis function
go_enrich_net(
  go_anno,
  go_deg_fc,
  padjust_method = "fdr",
  pvalue_cutoff = 0.5,
  qvalue_cutoff = 0.5,
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
data(kegg_anno)
head(kegg_anno)

data(kegg_deg_fc)
head(kegg_deg_fc)

# 2. Run go_enrich analysis function
res <- kegg_enrich(
  kegg_anno,
  kegg_deg_fc,
  padjust_method = "fdr",
  pvalue_cutoff = 1,
  qvalue_cutoff = 1
)
head(res)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::kegg_enrich

## ----kegg_enrich_bar, warning=FALSE, fig.width=10.00, fig.height=6.18---------
# 1. Load example datasets
data(kegg_anno)
# head(kegg_anno)

data(kegg_deg_fc)
# head(kegg_deg_fc)

# 2. Run kegg_enrich_bar analysis function
kegg_enrich_bar(
  kegg_anno,
  kegg_deg_fc,
  padjust_method = "fdr",
  pvalue_cutoff = 1,
  qvalue_cutoff = 1,
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
data(kegg_anno)
# head(kegg_anno)

data(kegg_deg_fc)
# head(kegg_deg_fc)

# 2. Run kegg_enrich_dot analysis function
kegg_enrich_dot(
  kegg_anno,
  kegg_deg_fc,
  padjust_method = "fdr",
  pvalue_cutoff = 1,
  qvalue_cutoff = 1,
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

## ----kegg_enrich_tree, warning=FALSE, fig.width=10.00, fig.height=6.18, eval=FALSE----
#  # 1. Load example datasets
#  data(kegg_anno)
#  # head(kegg_anno)
#  
#  data(kegg_deg_fc)
#  # head(kegg_deg_fc)
#  
#  # 2. Run kegg_enrich_tree analysis function
#  kegg_enrich_tree(
#    kegg_anno,
#    kegg_deg_fc,
#    padjust_method = "fdr",
#    pvalue_cutoff = 1,
#    qvalue_cutoff = 1,
#    sign_by = "p.adjust",
#    category_num = 20,
#    font_size = 4,
#    low_color = "#ff0000aa",
#    high_color = "#008800aa",
#    hclust_method = "complete",
#    ggTheme = "theme_light"
#  )

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::kegg_enrich_tree

## ----kegg_enrich_net, warning=FALSE, fig.width=10.00, fig.height=6.18---------
# 1. Load example datasets
data(kegg_anno)
# head(kegg_anno)

data(kegg_deg_fc)
# head(kegg_deg_fc)

# 2. Run kegg_enrich_net analysis function
kegg_enrich_net(
  kegg_anno,
  kegg_deg_fc,
  padjust_method = "fdr",
  pvalue_cutoff = 1,
  qvalue_cutoff = 1,
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
# 1. Load table_split_data example datasets
data(table_split_data)
head(table_split_data)

# 2. Run table_split plot function
res <- table_split(table_split_data, 
                  grouped_var = "variable", 
                  miss_drop = TRUE
                  )
head(res)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::table_split

## ----table_merge, warning=FALSE, fig.width=10.00, fig.height=6.18-------------
# 1. Load example datasets
data(table_merge_data)
head(table_merge_data)

# 2. Run function
res <- table_merge(
  table_merge_data,
  merge_vars = c("Ozone", "Solar.R", "Wind", "Temp"),
  new_var = "Variable",
  new_value = "Value",
  na_remove = FALSE
)
head(res)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::table_merge

## ----table_filter, warning=FALSE, fig.width=10.00, fig.height=6.18------------
# 1. Load example datasets
data(table_filter_data)
head(table_filter_data)

# 2. Run function
res <- table_filter(table_filter_data, 
                    height > 100 & eye_color == "black"
                    )
head(res)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::table_filter

## ----table_cross, warning=FALSE, fig.width=10.00, fig.height=6.18-------------
# 1. Load example datasets
data(table_cross_data1)
head(table_cross_data1)

data(table_cross_data2)
head(table_cross_data2)

# 2. Run function
res <- table_cross(
  table_cross_data1,
  table_cross_data2,
  inter_var = "geneID",
  left_index = TRUE,
  right_index = FALSE
)
head(res)

## -----------------------------------------------------------------------------
# Get help with command in R console.
# ?TOmicsVis::table_cross

