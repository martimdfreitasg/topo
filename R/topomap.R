#' Plot Alignment Matrix and Clusters
#'
#' @param fasta_path Path to the alignment fasta file or an XStringSet object.
#' @param metric Distance metric ("gaps", "residues", "both").
#' @param highlight_regions List of numeric vectors indicating column indices to highlight.
#' @param highlight_colors Character vector of hex codes corresponding to highlight_regions.
#' @param seed Random seed for reproducible clustering.
#' @param dist_mat Optional: Pre-calculated distance matrix from topodist().
#' @param cluster_mapping Optional: Pre-calculated mapping from topoclust().
#' @return A ggtree plot object containing the alignment heatmap and cluster ribbon.
#' @export
topomap <- function(fasta_path, metric = "gaps", highlight_regions = NULL, highlight_colors = NULL, seed = 12345, dist_mat = NULL, cluster_mapping = NULL) {

  # Retrieve upstream data prerequisites
  if (is.null(dist_mat)) dist_mat <- topodist(fasta_path, metric = metric)
  if (is.null(cluster_mapping)) cluster_mapping <- topoclust(fasta_path, metric = metric, seed = seed, dist_mat = dist_mat)

  unique_bin <- attr(dist_mat, "matrix")
  aln_dend <- stats::hclust(dist_mat, method = "average")

  tree <- ape::as.phylo(aln_dend)
  tree <- ape::reorder.phylo(tree, "cladewise")

  # Isolate unique cluster assignments for the ribbon mapping
  temp_anno <- cluster_mapping[!duplicated(cluster_mapping$Assigned_Group), ]
  cluster_anno <- data.frame(Cluster = temp_anno$Cluster)
  rownames(cluster_anno) <- temp_anno$Assigned_Group

  # Define base palette (0 = Residue, 1 = Gap)
  base_colors <- c("0" = "#333333", "1" = "#ffffff")

  # Overlay structural domain highlights onto the binary matrix
  if (!is.null(highlight_regions)) {
    if (is.null(highlight_colors) || length(highlight_colors) < length(highlight_regions)) {
      highlight_colors <- grDevices::hcl.colors(length(highlight_regions), palette = "Set 2")
    }
    for (i in seq_along(highlight_regions)) {
      cols <- highlight_regions[[i]]
      val <- i + 1
      temp_region <- unique_bin[, cols, drop = FALSE]
      temp_region[temp_region == 0] <- val
      unique_bin[, cols] <- temp_region
      base_colors[as.character(val)] <- highlight_colors[i]
    }
  }

  # Coerce matrix to dataframe for gheatmap compatibility
  aln_df <- as.data.frame(unique_bin)
  aln_df[] <- lapply(aln_df, as.character)

  # Generate discrete color palette for cluster assignments
  unique_clusters <- unique(as.character(cluster_anno$Cluster))
  cluster_colors <- grDevices::hcl.colors(length(unique_clusters), palette = "Pastel 1")
  names(cluster_colors) <- unique_clusters

  # Calculate dynamic structural offsets relative to the tree width
  tree_width <- max(ape::node.depth.edgelength(tree))
  clust_offset <- tree_width * 0.025
  mat_offset <- tree_width * 0.20
  text_offset <- tree_width * 3.25

  # Calculate dynamic right margin based on the longest sequence name
  max_label_length <- max(nchar(tree$tip.label))
  dynamic_right_margin <- (max_label_length * 5.5) + 15

  # Initialize invisible base tree geometry
  p_blank_tree <- ggtree::ggtree(tree, color = NA) +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))

  p_final <- suppressMessages({

    # Attach cluster ribbon layer
    p_clust <- ggtree::gheatmap(p_blank_tree, cluster_anno, offset = clust_offset, width = 0.1, colnames = FALSE, color = NA) +
      ggplot2::scale_fill_manual(values = cluster_colors, name = "Cluster")

    # Attach primary sequence alignment heatmap
    p_matrix <- ggtree::gheatmap(p_clust + ggnewscale::new_scale_fill(), aln_df, offset = mat_offset, width = 3, colnames = FALSE, color = NA) +
      ggplot2::scale_fill_manual(values = base_colors, guide = "none")

    # Finalize coordinate spacing and legend architecture
    p_matrix +
      ggtree::geom_tiplab(align = TRUE, offset = text_offset, size = 3.5, linetype = "blank") +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::theme(plot.margin = ggplot2::margin(10, dynamic_right_margin, 40, 10), # Implements the auto-sizing margin
                     legend.position = c(0.70, 0.0),
                     legend.justification = c(0.5, 0),
                     legend.direction = "horizontal",
                     legend.background = ggplot2::element_blank(),
                     legend.title = ggplot2::element_text(face = "bold"))
  })

  return(p_final)
}
