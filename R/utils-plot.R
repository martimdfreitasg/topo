#' Internal Helper to Build Cluster Ribbon and Alignment Plot
#'
#' @param p_base Base `ggtree` plot object.
#' @param tree `ape::phylo` object matching `p_base`.
#' @param unique_bin Deduplicated binary alignment matrix (`0` = residue, `1` = gap).
#' @param cluster_mapping Mapping `data.frame` from `topoclust()`.
#' @param highlight_regions Optional list of integer column index vectors.
#' @param highlight_colors Optional character vector of colors.
#' @return A `ggtree`/`ggplot` object.
#' @noRd
.build_topo_heatmap = function(p_base,
                                tree,
                                unique_bin,
                                cluster_mapping,
                                highlight_regions = NULL,
                                highlight_colors = NULL) {

  # Isolate unique clusters for mapping
  temp_anno = cluster_mapping[!duplicated(cluster_mapping$Assigned_Group), , drop = FALSE]
  cluster_anno = data.frame(
    Cluster = as.character(temp_anno$Cluster),
    stringsAsFactors = FALSE
  )
  rownames(cluster_anno) = temp_anno$Assigned_Group

  # Define base palette (0 = Residue, 1 = Gap)
  base_colors = c("0" = "#333333", "1" = "#ffffff")
  plot_bin = unique_bin
  n_cols = ncol(plot_bin)

  # Validate and overlay structural domain highlights onto the alignment
  if (!is.null(highlight_regions)) {
    if (!is.list(highlight_regions)) {
      stop("'highlight_regions' must be a list of positive integer column index vectors.", call. = FALSE)
    }

    if (is.null(highlight_colors)) {
      highlight_colors = grDevices::hcl.colors(length(highlight_regions), palette = "Set 2")
    } else if (length(highlight_colors) < length(highlight_regions)) {
      warning(
        "Fewer 'highlight_colors' supplied than 'highlight_regions'; falling back to default palette.",
        call. = FALSE
      )
      highlight_colors = grDevices::hcl.colors(length(highlight_regions), palette = "Set 2")
    }

    highlighted_cols = integer(0)
    for (i in seq_along(highlight_regions)) {
      cols = highlight_regions[[i]]
      if (!is.numeric(cols) || length(cols) == 0L || anyNA(cols) || any(cols != as.integer(cols))) {
        stop(
          sprintf("highlight_regions[[%d]] must be a non-empty integer vector of column indices.", i),
          call. = FALSE
        )
      }
      cols = as.integer(cols)
      if (any(cols < 1L | cols > n_cols)) {
        stop(
          sprintf(
            "highlight_regions[[%d]] contains out-of-bounds column indices (alignment has %d columns).",
            i, n_cols
          ),
          call. = FALSE
        )
      }

      overlap = intersect(highlighted_cols, cols)
      if (length(overlap) > 0L) {
        warning(
          sprintf(
            "highlight_regions[[%d]] overlaps with an earlier region at %d position(s); the earlier region takes precedence.",
            i, length(overlap)
          ),
          call. = FALSE
        )
      }
      highlighted_cols = union(highlighted_cols, cols)

      val = i + 1L
      temp_region = plot_bin[, cols, drop = FALSE]
      temp_region[temp_region == 0] = val
      plot_bin[, cols] = temp_region
      base_colors[as.character(val)] = highlight_colors[i]
    }
  }

  # Coerce matrix to character dataframe for gheatmap compatibility
  aln_df = as.data.frame(plot_bin, stringsAsFactors = FALSE)
  aln_df[] = lapply(aln_df, as.character)

  # Discrete cluster palette: numeric order for clusters, neutral gray for "Unassigned"
  unique_clusters = unique(cluster_anno$Cluster)
  has_unassigned = "Unassigned" %in% unique_clusters
  assigned_clusters = setdiff(unique_clusters, "Unassigned")

  if (length(assigned_clusters) > 0L) {
    clust_nums = suppressWarnings(as.numeric(sub("^Cluster\\s+", "", assigned_clusters)))
    if (!anyNA(clust_nums)) {
      assigned_clusters = assigned_clusters[order(clust_nums)]
    } else {
      assigned_clusters = sort(assigned_clusters)
    }
    cluster_colors = grDevices::hcl.colors(length(assigned_clusters), palette = "Pastel 1")
    names(cluster_colors) = assigned_clusters
  } else {
    cluster_colors = character(0)
  }

  if (has_unassigned) {
    cluster_colors = c(cluster_colors, Unassigned = "#BDBDBD")
  }

  cluster_anno$Cluster = factor(cluster_anno$Cluster, levels = names(cluster_colors))

  # Dynamic offsets relative to tree width
  tree_width = max(ape::node.depth.edgelength(tree))
  clust_offset = tree_width * 0.025
  mat_offset = tree_width * 0.20
  text_offset = tree_width * 3.25

  # Dynamic right margin based on longest sequence label
  max_label_length = max(nchar(tree$tip.label))
  dynamic_right_margin = (max_label_length * 5.5) + 15

  p_final = suppressMessages({
    p_clust = ggtree::gheatmap(
      p_base,
      cluster_anno,
      offset = clust_offset,
      width = 0.1,
      colnames = FALSE,
      color = NA
    ) +
      ggplot2::scale_fill_manual(values = cluster_colors, name = "Cluster")

    p_matrix = ggtree::gheatmap(
      p_clust + ggnewscale::new_scale_fill(),
      aln_df,
      offset = mat_offset,
      width = 3,
      colnames = FALSE,
      color = NA
    ) +
      ggplot2::scale_fill_manual(values = base_colors, guide = "none")

    p_matrix +
      ggtree::geom_tiplab(align = TRUE, offset = text_offset, size = 3.5, linetype = "blank") +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::theme(
        plot.margin = ggplot2::margin(10, dynamic_right_margin, 40, 10),
        legend.position = "inside",
        legend.position.inside = c(0.70, 0.0),
        legend.justification = c(0.5, 0),
        legend.direction = "horizontal",
        legend.background = ggplot2::element_blank(),
        legend.title = ggplot2::element_text(face = "bold")
      )
  })

  return(p_final)
}
