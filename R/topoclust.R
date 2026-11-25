#' Extract Dynamic Topological Cluster Assignments
#'
#' Performs average-linkage (UPGMA) hierarchical clustering on alignment
#' topological distances and identifies discrete clusters via dynamic tree
#' cutting.
#'
#' @details
#' Rather than applying a single static height threshold, which performs poorly
#' when sequence clades diverge unevenly or have unbalanced sizes, clusters are
#' made from dendrogram branch shape using
#' `dynamicTreeCut::cutreeDynamic()` with the hybrid method (`pamStage = TRUE`).
#'
#' Sequences labeled `"Unassigned"` represent topologies that could not be
#' grouped into a cluster meeting `minClusterSize` distinct topologies; this
#' indicates an isolated structural pattern rather than a data error.
#'
#' To obtain reproducible results if clustering ties occur, call
#' `set.seed()` prior to running `topoclust()`.
#'
#' @param fasta_path Path to a FASTA alignment file or an `ape::AAbin` /
#'   `ape::DNAbin` object. Optional if `dist_mat` is provided.
#' @param metric Distance metric passed to `topodist()` when `dist_mat` is
#'   `NULL`: one of `"gaps"` (default), `"residues"`, or `"both"`.
#' @param dist_mat Optional pre-calculated `dist` object from `topodist()`.
#' @param deepSplit Controls cluster splitting sensitivity in
#'   `dynamicTreeCut::cutreeDynamic()`. Accepts an integer in `0:4` (higher
#'   values produce more, smaller clusters) or a logical. Default is `2`.
#' @param minClusterSize Minimum number of distinct sequence topologies
#'   required to form a cluster. Default is `2`.
#'
#' @return A `data.frame` with three columns:
#'   \itemize{
#'     \item `Original_Sequence`: Original sequence header from the input.
#'     \item `Assigned_Group`: Deduplicated topology group or singleton name.
#'     \item `Cluster`: Factor indicating cluster membership (`"Cluster 1"`,
#'       `"Cluster 2"`, ..., and `"Unassigned"` if present).
#'   }
#'
#' @references
#' Langfelder, P., Zhang, B., & Horvath, S. (2008). Defining clusters from a
#' hierarchical cluster tree: the Dynamic Tree Cut package for R.
#' *Bioinformatics*, 24(5), 719--720. \doi{10.1093/bioinformatics/btm563}
#'
#' @export
#' @examples
#' fasta_path = system.file("extdata", "test.fasta", package = "topo")
#' set.seed(12345)
#' clusters = topoclust(fasta_path, metric = "gaps")
#' head(clusters)
topoclust = function(fasta_path = NULL,
                      metric = c("gaps", "residues", "both"),
                      dist_mat = NULL,
                      deepSplit = 2,
                      minClusterSize = 2) {

  metric = match.arg(metric)

  if (!(is.logical(deepSplit) && length(deepSplit) == 1L && !is.na(deepSplit)) &&
      !(is.numeric(deepSplit) && length(deepSplit) == 1L && !is.na(deepSplit) &&
        deepSplit >= 0 && deepSplit <= 4)) {
    stop("'deepSplit' must be a single logical value or a number between 0 and 4.", call. = FALSE)
  }
  if (!is.numeric(minClusterSize) || length(minClusterSize) != 1L || is.na(minClusterSize) ||
      minClusterSize < 1 || minClusterSize != as.integer(minClusterSize)) {
    stop("'minClusterSize' must be a single positive integer.", call. = FALSE)
  }

  if (is.null(dist_mat)) {
    dist_mat = topodist(fasta_path, metric = metric)
  }
  if (!inherits(dist_mat, "dist")) {
    stop("'dist_mat' must be a 'dist' object produced by topodist().", call. = FALSE)
  }
  if (attr(dist_mat, "Size") < 2L) {
    stop("At least 2 distinct topologies are required for clustering.", call. = FALSE)
  }

  unique_bin = attr(dist_mat, "matrix")
  row_ids = if (!is.null(unique_bin)) rownames(unique_bin) else labels(dist_mat)

  aln_dend = stats::hclust(dist_mat, method = "average")

  dynamic_mods = dynamicTreeCut::cutreeDynamic(
    dendro = aln_dend,
    distM = as.matrix(dist_mat),
    deepSplit = deepSplit,
    pamStage = TRUE,
    minClusterSize = as.integer(minClusterSize),
    verbose = 0
  )

  cluster_labels = ifelse(dynamic_mods == 0L, "Unassigned", paste("Cluster", dynamic_mods))
  assigned_nums = sort(unique(dynamic_mods[dynamic_mods > 0L]))
  cluster_levels = c(
    if (length(assigned_nums) > 0L) paste("Cluster", assigned_nums),
    if (any(dynamic_mods == 0L)) "Unassigned"
  )

  cluster_vec = stats::setNames(
    factor(cluster_labels, levels = cluster_levels),
    row_ids
  )

  assigned_groups = attr(dist_mat, "assigned_groups")
  if (is.null(assigned_groups)) assigned_groups = row_ids
  original_names = attr(dist_mat, "original_names")
  if (is.null(original_names)) original_names = row_ids

  mapping_df = data.frame(
    Original_Sequence = unname(original_names),
    Assigned_Group = unname(assigned_groups),
    Cluster = unname(cluster_vec[assigned_groups]),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  return(mapping_df)
}
