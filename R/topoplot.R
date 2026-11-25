#' Plot Combined Topological Dendrogram, Clusters, and Alignment
#'
#' Wrapper that computes alignment topological distances, dynamic
#' tree-cut clusters, and an optionally bootstrapped UPGMA dendrogram, then
#' renders them together with the alignment and optional domain
#' highlights.
#'
#' @details
#' Runs `topodist()`, `topoclust()`, and `topodend()` once and
#' passes the cached objects through the plotting pipeline to avoid redundant
#' I/O or distance computations. See `topodist()` for metric definitions,
#' `topoclust()` for dynamic clustering parameters, and `topodend()` for
#' bootstrap and polytomy collapsing.
#'
#' For reproducible bootstrap values and clustering, call `set.seed()` prior to
#' running `topoplot()`.
#'
#' @param fasta_path Path to a FASTA alignment file or an `ape::AAbin` /
#'   `ape::DNAbin` object.
#' @param metric Distance metric (`"gaps"`, `"residues"`, or `"both"`).
#' @param highlight_regions Optional list of positive integer vectors defining
#'   alignment column indices to highlight.
#' @param highlight_colors Optional character vector of colors corresponding to
#'   `highlight_regions`. Defaults to `"Set 2"` from `grDevices::hcl.colors()`.
#' @param bootstrap Positive integer number of bootstrap replicates. When `1`
#'   (default), bootstrapping is omitted.
#' @param boot_threshold Optional numeric percentage in `[0, 100]` below which
#'   internal nodes are collapsed into polytomies.
#' @param deepSplit Cluster splitting sensitivity passed to `topoclust()`.
#'   Default is `2`.
#' @param minClusterSize Minimum cluster size passed to `topoclust()`. Default
#'   is `2`.
#'
#' @return A combined `ggtree`/`ggplot` object.
#'
#' @references
#' Langfelder, P., Zhang, B., & Horvath, S. (2008). Defining clusters from a
#' hierarchical cluster tree: the Dynamic Tree Cut package for R.
#' *Bioinformatics*, 24(5), 719--720. \doi{10.1093/bioinformatics/btm563}
#'
#' Yu, G., Smith, D. K., Zhu, H., Guan, Y., & Lam, T. T.-Y. (2017). ggtree: an
#' R package for visualization and annotation of phylogenetic trees with their
#' covariates and other associated data. *Methods in Ecology and Evolution*,
#' 8(1), 28--36. \doi{10.1111/2041-210X.12628}
#'
#' @export
#' @examples
#' fasta_path = system.file("extdata", "test.fasta", package = "topo")
#' final_plot = topoplot(fasta_path, metric = "gaps")
#' final_plot
topoplot = function(fasta_path = NULL,
                     metric = c("gaps", "residues", "both"),
                     highlight_regions = NULL,
                     highlight_colors = NULL,
                     bootstrap = 1,
                     boot_threshold = NULL,
                     deepSplit = 2,
                     minClusterSize = 2) {

  metric = match.arg(metric)

  dist_mat = topodist(fasta_path, metric = metric)
  cluster_mapping = topoclust(
    fasta_path = fasta_path,
    metric = metric,
    dist_mat = dist_mat,
    deepSplit = deepSplit,
    minClusterSize = minClusterSize
  )

  dend_obj = topodend(
    fasta_path = fasta_path,
    metric = metric,
    dist_mat = dist_mat,
    bootstrap = bootstrap,
    show_tiplab = FALSE,
    boot_threshold = boot_threshold
  )

  .build_topo_heatmap(
    p_base = dend_obj$as_ggtree,
    tree = dend_obj$as_phylo,
    unique_bin = attr(dist_mat, "matrix"),
    cluster_mapping = cluster_mapping,
    highlight_regions = highlight_regions,
    highlight_colors = highlight_colors
  )
}
