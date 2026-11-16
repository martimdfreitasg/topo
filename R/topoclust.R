#' Extract Cluster and Group Mapping
#'
#' @param fasta_path Path to the alignment fasta file or an XStringSet object.
#' @param metric Distance metric ("gaps", "residues", "both").
#' @param dist_mat Optional: Pre-calculated distance matrix from topodist().
#' @return A dataframe mapping original sequences to their deduplicated group and assigned cluster.
#' @export
#' @examples
#' fasta_path <- system.file("extdata", "test.fasta", package = "topo")
#' clusters <- topoclust(fasta_path, metric = "gaps")
topoclust <- function(fasta_path, metric = "gaps", dist_mat = NULL) {

  # Retrieve or compute the distance matrix
  if (is.null(dist_mat)) {
    dist_mat <- topodist(fasta_path, metric = metric)
  }
  unique_bin <- attr(dist_mat, "matrix")

  # Perform Average Linkage (UPGMA) hierarchical clustering
  aln_dend <- stats::hclust(dist_mat, method = "average")

  # Define discrete clusters using dynamic tree cutting
  dynamic_mods <- dynamicTreeCut::cutreeDynamic(
    dendro = aln_dend,
    distM = as.matrix(dist_mat),
    deepSplit = 2,
    pamStage = TRUE,
    minClusterSize = 2,
    verbose = 0
  )

  # Format cluster assignments and append to the unique sequence data
  cluster_labels <- ifelse(dynamic_mods == 0, "Unassigned", paste("Cluster", dynamic_mods))
  cluster_anno <- data.frame(Cluster = as.factor(cluster_labels))
  rownames(cluster_anno) <- rownames(unique_bin)

  # Retrieve sequence mapping directly from the distance matrix attributes
  assigned_groups <- attr(dist_mat, "assigned_groups")
  original_names <- attr(dist_mat, "original_names")

  # Assemble the final mapping frame linking original sequences to their clusters
  mapping_df <- data.frame(
    Original_Sequence = original_names,
    Assigned_Group = assigned_groups,
    Cluster = cluster_anno[assigned_groups, "Cluster"],
    row.names = NULL
  )

  return(mapping_df)
}
