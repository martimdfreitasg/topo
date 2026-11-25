#' Plot Topological Distance Dendrogram
#'
#' Constructs an average-linkage (UPGMA) dendrogram from alignment topological
#' distances, with optional column-resampling bootstrap support and
#' polytomy condensation for weakly supported nodes.
#'
#' @details
#' When `bootstrap > 1`, alignment columns are resampled with replacement via
#' `ape::boot.phylo()`. Because `ape::boot.phylo()` returns raw replicate
#' counts, values are automatically converted to percentages (`0` to `100`)
#' before labeling internal nodes and comparing against `boot_threshold`.
#'
#' When `boot_threshold` is set, internal edges leading to nodes with bootstrap
#' support strictly below `boot_threshold` are collapsed to length `0` in
#' root-to-tip (`"cladewise"`) order, and each collapsed edge's length is
#' redistributed to its immediate child edges. This preserves the ultrametric
#' tip depth of the UPGMA dendrogram so the distance scale axis remains exact
#' for all tips and surviving nodes.
#'
#' For reproducible bootstrap values, call `set.seed()` prior to running
#' `topodend()`.
#'
#' @param fasta_path Path to a FASTA alignment file or an `ape` `AAbin`/`DNAbin`
#' object.
#'
#'   Optional if `dist_mat` or `topomatrix` is provided.
#' @param topomatrix Optional numeric binary matrix passed to `topodist()`.
#' @param metric Distance metric (`"gaps"`, `"residues"`, or `"both"`). Ignored
#'   when `dist_mat` carries a `"metric"` attribute from `topodist()`.
#' @param dist_mat Optional pre-calculated `dist` object from `topodist()`.
#' @param bootstrap Positive integer number of bootstrap replicates. When `1`
#'   (default), bootstrapping is omitted.
#' @param show_tiplab Logical indicating whether to display tip labels on the
#'   dendrogram. Default is `TRUE`.
#' @param boot_threshold Optional numeric percentage in `[0, 100]`. Internal
#'   nodes with bootstrap support strictly below this threshold are collapsed
#'   into polytomies.
#'
#' @return A list with two elements:
#'   \itemize{
#'     \item `as_phylo`: The `ape::phylo` dendrogram (with bootstrap percentage
#'       `node.label`s and collapsed edges if requested).
#'     \item `as_ggtree`: The rendered `ggtree`/`ggplot` object.
#'   }
#'
#' @references
#' Felsenstein, J. (1985). Confidence limits on phylogenies: an approach using
#' the bootstrap. *Evolution*, 39(4), 783--791.
#' \doi{10.1111/j.1558-5646.1985.tb00420.x}
#'
#' Yu, G., Smith, D. K., Zhu, H., Guan, Y., & Lam, T. T.-Y. (2017). ggtree: an
#' R package for visualization and annotation of phylogenetic trees with their
#' covariates and other associated data. *Methods in Ecology and Evolution*,
#' 8(1), 28--36. \doi{10.1111/2041-210X.12628}
#'
#' @export
#' @examples
#' fasta_path = system.file("extdata", "test.fasta", package = "topo")
#' dend = topodend(fasta_path, metric = "gaps")
#' dend$as_ggtree
topodend = function(fasta_path = NULL,
                     topomatrix = NULL,
                     metric = c("gaps", "residues", "both"),
                     dist_mat = NULL,
                     bootstrap = 1,
                     show_tiplab = TRUE,
                     boot_threshold = NULL) {

  metric = match.arg(metric)

  if (!is.numeric(bootstrap) || length(bootstrap) != 1L || is.na(bootstrap) ||
      bootstrap < 1 || bootstrap != as.integer(bootstrap)) {
    stop("'bootstrap' must be a single positive integer.", call. = FALSE)
  }
  if (!is.null(boot_threshold)) {
    if (!is.numeric(boot_threshold) || length(boot_threshold) != 1L || is.na(boot_threshold) ||
        boot_threshold < 0 || boot_threshold > 100) {
      stop("'boot_threshold' must be a single numeric percentage between 0 and 100.", call. = FALSE)
    }
  }

  if (is.null(dist_mat)) {
    dist_mat = topodist(fasta_path, topomatrix = topomatrix, metric = metric)
  }
  if (!inherits(dist_mat, "dist")) {
    stop("'dist_mat' must be a 'dist' object produced by topodist().", call. = FALSE)
  }
  if (attr(dist_mat, "Size") < 2L) {
    stop("At least 2 distinct topologies are required to build a dendrogram.", call. = FALSE)
  }

  metric_used = attr(dist_mat, "metric")
  if (is.null(metric_used)) {
    metric_used = metric
  }

  aln_dend = stats::hclust(dist_mat, method = "average")
  tree = ape::as.phylo(aln_dend)
  tree = ape::reorder.phylo(tree, "cladewise")
  Ntip = length(tree$tip.label)

  # Optional bootstrapping and polytomy condensation
  if (bootstrap > 1L) {
    bin_mat = attr(dist_mat, "matrix")
    if (is.null(bin_mat)) {
      stop(
        "Bootstrapping requires the binary matrix attribute attached to 'dist_mat' by topodist().",
        call. = FALSE
      )
    }

    phylo_func = function(x) {
      dist_mat_temp = topodist(fasta_path = NULL, topomatrix = x, metric = metric_used)
      aln_dend_temp = stats::hclust(dist_mat_temp, method = "average")
      ape::as.phylo(aln_dend_temp)
    }

    boot_counts = ape::boot.phylo(
      phy = tree,
      x = bin_mat,
      FUN = phylo_func,
      B = as.integer(bootstrap),
      quiet = TRUE
    )

    boot_pct = (boot_counts / bootstrap) * 100
    tree$node.label = ifelse(is.na(boot_pct), "", as.character(round(boot_pct)))

    if (!is.null(boot_threshold)) {
      failed_nodes = which(!is.na(boot_pct) & boot_pct < boot_threshold)

      if (length(failed_nodes) > 0L) {
        nodes_to_collapse = Ntip + failed_nodes

        # Walk edges in preorder ("cladewise") so parent collapses propagate
        # cleanly into child edges, keeping all tips at tree_width
        for (e in seq_len(nrow(tree$edge))) {
          child_node = tree$edge[e, 2L]
          if (child_node %in% nodes_to_collapse) {
            collapsed_len = tree$edge.length[e]
            tree$edge.length[e] = 0
            child_edges = which(tree$edge[, 1L] == child_node)
            tree$edge.length[child_edges] = tree$edge.length[child_edges] + collapsed_len
          }
        }

        tree$node.label[failed_nodes] = ""
      }
    }
  }

  # Unified dendrogram and distance scale
  y_axis = 0
  y_ticks = -(Ntip * 0.02)
  y_labels = -(Ntip * 0.05)
  y_title = -(Ntip * 0.10)
  tree_width = max(ape::node.depth.edgelength(tree))
  max_dist = 2 * tree_width
  tick_step = ifelse(metric_used == "both", 0.1, 0.2)
  dist_ticks = seq(0, max_dist, by = tick_step)
  if (length(dist_ticks) < 2L && max_dist > 0) {
    dist_ticks = c(0, max_dist)
  }
  axis_breaks = tree_width - (dist_ticks / 2)
  axis_labels = sprintf("%.2f", dist_ticks)

  axis_title = switch(
    metric_used,
    gaps = "Distance (Jaccard - Gaps)",
    residues = "Distance (Jaccard - Residues)",
    both = "Distance (Hamming)"
  )

  p_tree = suppressMessages({
    ggtree::ggtree(tree, linewidth = 0.5) +
      ggtree::theme_tree() +
      ggplot2::annotate("segment", x = 0, xend = tree_width, y = y_axis, yend = y_axis, linewidth = 0.4) +
      ggplot2::annotate("segment", x = axis_breaks, xend = axis_breaks, y = y_axis, yend = y_ticks, linewidth = 0.4) +
      ggplot2::annotate("text", x = axis_breaks, y = y_labels, label = axis_labels, size = 3) +
      ggplot2::annotate("text", x = tree_width / 2, y = y_title, label = axis_title, size = 3.5) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.15, 0.05)))
  })

  if (bootstrap > 1L) {
    p_tree = p_tree +
      ggplot2::geom_label(
        ggplot2::aes(label = label),
        data = function(d) d[!d$isTip & !is.na(d$label) & d$label != "", ],
        size = 3
      )
  }

  if (show_tiplab) {
    p_tree = p_tree + ggtree::geom_tiplab(size = 3)
  }

  return(list(as_phylo = tree, as_ggtree = p_tree))
}
