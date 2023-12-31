#' @title KEGG enrichment analysis and net plot (None/Exist Reference Genome).
#' @description KEGG enrichment analysis and net plot (None/Exist Reference Genome).
#' @author benben-miao
#'
#' @return Plot: KEGG enrichment analysis and net plot (None/Exist Reference Genome).
#' @param kegg_anno Dataframe: GO and KEGG annotation of background genes (1st-col: Genes, 2nd-col: biological_process, 3rd-col: cellular_component, 4th-col: molecular_function, 5th-col: kegg_pathway).
#' @param degs_list Dataframe: degs list.
#' @param padjust_method Character: P-value adjust to Q-value. Default: "fdr" (false discovery rate), options: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' @param pvalue_cutoff Numeric: P-value cutoff. Recommend: small than 0.05.
#' @param qvalue_cutoff Numeric: Q-value cutoff. Recommend: small than 0.05.
#' @param category_num Numeric: categories number to display. Default: 20, min: 1, max: NULL.
#' @param net_layout Character: network layout. Default: "circle", options: 'star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 'randomly', 'fr', 'kk', 'drl' or 'lgl'.
#' @param net_circular Logical: network circular. Default: TRUE, options: TRUE, FALSE.
#' @param low_color Character: low value (p-value or q-value) color (color name or hex value).
#' @param high_color Character: high value (p-value or q-value) color (color name or hex value).
#'
#' @import ggplot2
#' @import ggsci
#' @importFrom reshape2 melt
#' @importFrom tidyr separate_rows separate drop_na
#' @importFrom clusterProfiler enricher
#' @importFrom dplyr distinct
#' @import enrichplot
#' @import ggnewscale
#' @export
#'
#' @examples
#' # 1. Library TOmicsVis package
#' library(TOmicsVis)
#'
#' # 2. Use example dataset
#' data(gene_go_kegg)
#' head(gene_go_kegg)
#'
#' # 3. Default parameters
#' kegg_enrich_net(gene_go_kegg[,c(1,5)], gene_go_kegg[100:200,1])
#'
#' # 4. Set category_num = 10
#' kegg_enrich_net(gene_go_kegg[,c(1,5)], gene_go_kegg[100:200,1], category_num = 10)
#'
kegg_enrich_net <- function(kegg_anno,
														degs_list,
													 padjust_method = "fdr",
													 pvalue_cutoff = 0.05,
													 qvalue_cutoff = 0.05,
													 category_num = 20,
													 net_layout = "circle",
													 net_circular = TRUE,
													 low_color = "#ff0000aa",
													 high_color = "#008800aa"
													){
	# -> 2. Data Parameters
	# padjust_method <- "fdr"
	# ChoiceBox: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"

	# pvalue_cutoff <- 0.30
	# Slider: 0.30, 0.00, 0.01, 1.00

	# qvalue_cutoff <- 0.50
	# Slider: 0.50, 0.00, 0.01, 1.00
	# <- 2. Data Parameters

	# -> 3. Data
	gene_kegg <- kegg_anno
	degs_list <- degs_list

	# deg_fc["log2FC"] <- 2^(deg_fc["log2FC"])
	# deg_list <- with(deg_fc, setNames(log2FC, id))

	gene_kegg7 <- separate_rows(data = gene_kegg,
															"kegg_pathway",
															sep = ";"
	)

	gene_kegg8 <- separate(gene_kegg7,
												 "kegg_pathway",
												 c("kegg_pathway", "description"),
												 "\\("
	)

	gene_kegg9 <- drop_na(gene_kegg8)
	gene_kegg9["description"] <- gsub(")", "", gene_kegg9$description)

	enrich_kegg <- enricher(gene = degs_list,
													TERM2GENE = data.frame(gene_kegg9[,2],gene_kegg9[,1]),
													TERM2NAME = gene_kegg9[,2:3],
													pvalueCutoff = pvalue_cutoff,
													pAdjustMethod = padjust_method, # "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
													qvalueCutoff = qvalue_cutoff,
													minGSSize = 1,
													maxGSSize = 1000
	)

	enrich_result <- enrich_kegg@result

	# write.table(enrich_result,
	# 						file = "Results.txt",
	# 						append = FALSE,
	# 						sep = "\t",
	# 						quote = TRUE,
	# 						na = "NA"
	# )
	# <- 3. Data

	# -> 4. Plot parameters
	# fonts <- "Times"
	# ChoiceBox: "Times", "Palatino", "Bookman", "Courier", "Helvetica", "URWGothic", "NimbusMon", "NimbusSan"

	# category_num <- 30
	# ChoiceBox: 30, 10, 1, 50

	# low_color <- "#ff0000aa"
	# ColorPicker

	# high_color <- "#0000ffaa"
	# ColorPicker

	# net_layout <- "circle"
	# ChoiceBox: 'star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 'randomly', 'fr', 'kk', 'drl' or 'lgl'.

	# isCircular <- "circular"
	# if (isCircular == "circular") {
	# 	circular <- TRUE
	# } else if (isCircular == "scattered") {
	# 	circular <- FALSE
	# }
	# # ChoiceBox: "circular", "scattered"


	cateLabelScale <- 0.8
	# Slider: 0.8, 0.1, 0.1, 2.0

	geneLabelScale <- 0.8
	# Slider: 0.8, 0.1, 0.1, 2.0
	# <- 4. Plot parameters

	# -> 5. Plot
	p <- cnetplot(
		x = enrich_kegg,
		showCategory = category_num,
		# foldChange = deg_list,
		layout = net_layout, # 'star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 'randomly', 'fr', 'kk', 'drl' or 'lgl'.
		# colorEdge = TRUE,
		circular = net_circular,
		node_label = "all",
		# cex_category = 1,
		# cex_gene = 1,
		# cex_label_category = cateLabelScale,
		# cex_label_gene = geneLabelScale,
		# color_category = "#ff0000",
		# color_gene = "#008000",
		shadowtext = "all",
		color.params = list(foldChange = NULL,
												edge = TRUE),
		cex.params = list(category_node = 1,
											gene_node = 1,
											category_label = cateLabelScale,
											gene_label = geneLabelScale)
	) +
		# guides(color = guide_legend(title="New Legend Title")) +
		labs(color = "Genes") +
		theme(
			# text = element_text(family = fonts)
		) +
		scale_fill_gradient(low = low_color, high = high_color,
												space = "Lab",
												guide = "colourbar", aesthetics = "fill")

	# p
	# <- 5. Plot

	return(p)
	invisible()
}
