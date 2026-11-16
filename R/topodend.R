#' Plot Distance Dendrogram
#'
#' @param fasta_path Path to the alignment fasta file or an XStringSet object.
#' @param topomatrix Optional: Pre-processed matrix to override fasta path.
#' @param metric Distance metric ("gaps", "residues", "both").
#' @param dist_mat Optional: Pre-calculated distance matrix from topodist().
#' @param bootstrap Optional: Bootstrap iterations. When set to 1, bootstrap is omitted.
#' @param seed Optional: Random seed for bootstrapping.
#' @param show_tiplab Optional: Logical to toggle tip labels on the base tree. Default is TRUE.
#' @param boot_threshold Optional: Numeric threshold below which nodes are collapsed into a polytomy.
#' @return A list containing the raw phylo object and the rendered ggtree object.
#' @export
topodend <- function(fasta_path, topomatrix = NULL, metric = "gaps", dist_mat = NULL, bootstrap = 1, seed = 1, show_tiplab = TRUE, boot_threshold = NULL) {

  if (is.null(dist_mat)) {
    dist_mat <- topodist(fasta_path, topomatrix, metric = metric)
  }

  metric_used <- attr(dist_mat, "metric")
  aln_dend <- stats::hclust(dist_mat, method = "average")

  # Convert to phylo format and reorder edges for downstream ggtree compatibility
  tree <- ape::as.phylo(aln_dend)
  tree <- ape::reorder.phylo(tree, "cladewise")

  # Calculate dynamic Y-axis spacing based on sequence volume
  Ntip <- length(tree$tip.label)
  y_axis <- 0
  y_ticks <- -(Ntip * 0.02)
  y_labels <- -(Ntip * 0.05)
  y_title <- -(Ntip * 0.10)

  # Process without bootstrapping
  if (bootstrap <= 1) {
    tree_width <- max(ape::node.depth.edgelength(tree))
    tick_step <- ifelse(metric_used == "both", 0.05, 0.1)
    axis_breaks <- seq(0, tree_width, by = tick_step)
    axis_labels <- sprintf("%.2f", abs(tree_width - axis_breaks))

    axis_title <- switch(metric_used,
                         "gaps" = "Distance (Jaccard - Gaps)",
                         "residues" = "Distance (Jaccard - Residues)",
                         "both" = "Distance (Hamming)")

    p_tree <- ggtree::ggtree(tree, linewidth = 0.5) +
      ggtree::theme_tree()

    # Append standard distance scale
    p_tree <- p_tree +
      ggplot2::annotate("segment", x = 0, xend = tree_width, y = y_axis, yend = y_axis, linewidth = 0.4) +
      ggplot2::annotate("segment", x = axis_breaks, xend = axis_breaks, y = y_axis, yend = y_ticks, linewidth = 0.4) +
      ggplot2::annotate("text", x = axis_breaks, y = y_labels, label = axis_labels, size = 3) +
      ggplot2::annotate("text", x = tree_width / 2, y = y_title, label = axis_title, size = 3.5) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.15, 0.05)))

    if (show_tiplab) {
      p_tree <- p_tree + ggtree::geom_tiplab(size = 3)
    }

    # Process with bootstrap iterations
  } else {
    phylo_func <- function(topomatrix) {
      dist_mat_temp <- topodist(fasta_path, topomatrix, metric = metric)
      aln_dend <- stats::hclust(dist_mat_temp, method = "average")
      return(ape::as.phylo(aln_dend))
    }

    boot_vals <- ape::boot.phylo(tree, attr(dist_mat, "matrix"), phylo_func, bootstrap)
    tree$node.label <- as.character(boot_vals)

    # Evaluate node support and generate visual polytomies for low-confidence branches
    if (!is.null(boot_threshold)) {
      failed_nodes <- which(boot_vals < boot_threshold)

      if (length(failed_nodes) > 0) {
        nodes_to_collapse <- Ntip + failed_nodes
        edges_to_collapse <- which(tree$edge[, 2] %in% nodes_to_collapse)

        # Collapse branches visually by setting length to zero
        tree$edge.length[edges_to_collapse] <- 0
        tree$node.label[failed_nodes] <- "" # Erase labels to prevent text overlap at the origin
      }
    }

    tree_width <- max(ape::node.depth.edgelength(tree))
    tick_step <- ifelse(metric == "both", 0.05, 0.1)
    axis_breaks <- seq(0, tree_width, by = tick_step)
    axis_labels <- sprintf("%.2f", abs(tree_width - axis_breaks))

    axis_title <- switch(metric,
                         "gaps" = "Distance (Jaccard - Gaps)",
                         "residues" = "Distance (Jaccard - Residues)",
                         "both" = "Distance (Hamming)")

    p_tree <- ggtree::ggtree(tree, linewidth = 0.5) +
      ggtree::theme_tree()

    # Append distance scale and map remaining bootstrap labels
    p_tree <- p_tree +
      ggplot2::annotate("segment", x = 0, xend = tree_width, y = y_axis, yend = y_axis, linewidth = 0.4) +
      ggplot2::annotate("segment", x = axis_breaks, xend = axis_breaks, y = y_axis, yend = y_ticks, linewidth = 0.4) +
      ggplot2::annotate("text", x = axis_breaks, y = y_labels, label = axis_labels, size = 3) +
      ggplot2::annotate("text", x = tree_width / 2, y = y_title, label = axis_title, size = 3.5) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.15, 0.05))) +
      ggplot2::geom_label(
        ggplot2::aes(label = label),
        data = function(d) d[!d$isTip & !is.na(d$label) & d$label != "", ], # Dynamic subsetting prevents missing column errors
        size = 3
      )

    if (show_tiplab) {
      p_tree <- p_tree + ggtree::geom_tiplab(size = 3)
    }
  }

  return(list(as_phylo = tree, as_ggtree = p_tree))
}
