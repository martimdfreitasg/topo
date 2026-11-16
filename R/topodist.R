#' Calculate Alignment Distances
#'
#' @param fasta_path Path to the alignment fasta file or an XStringSet object.
#' @param topomatrix Optional: Pre-processed matrix to override fasta file input.
#' @param metric Distance metric ("gaps", "residues", "both").
#' @return A dist object with the deduplicated matrix and metadata attached as attributes.
#' @export
#' @examples
#' fasta_path <- system.file("extdata", "test.fasta", package = "topo")
#' gap_dist <- topodist(fasta_path, metric = "gaps")
topodist <- function(fasta_path, topomatrix = NULL, metric = "gaps") {

  if (is.null(topomatrix)) {
    # Load sequence alignment from file path
    if (is.character(fasta_path)) {
      aln <- seqinr::read.fasta(fasta_path, seqtype = "AA", as.string = FALSE, whole.header = TRUE)
    } else if (inherits(fasta_path, "XStringSet")) {
      aln <- fasta_path
    } else {
      stop("Input must be a valid file path or a Biostrings XStringSet object.")
    }

    # Pad sequences with terminal gaps to ensure uniform alignment length for matrix conversion
    seq_lens <- lengths(aln)
    if (length(unique(seq_lens)) > 1) {
      max_len <- max(seq_lens)
      aln <- lapply(aln, function(x) c(x, rep("-", max_len - length(x))))
    }

    # Coerce alignment to a binary matrix (1 = Gap, 0 = Residue)
    bin_mat <- (do.call(rbind, aln) == "-") + 0

    # Deduplicate identical sequences into representative groups to optimize clustering
    row_sigs <- apply(bin_mat, 1, paste, collapse = "")
    unique_idx <- !duplicated(row_sigs)
    unique_bin <- bin_mat[unique_idx, , drop = FALSE]

    # Assign group names to sequences that share identical topologies
    is_group <- table(row_sigs)[row_sigs[unique_idx]] > 1
    rownames(unique_bin)[is_group] <- paste0("Group ", seq_len(sum(is_group)))

    # Store mapping dictionary to bypass redundant FASTA reading in downstream functions
    signature_to_group <- stats::setNames(rownames(unique_bin), row_sigs[unique_idx])
    assigned_groups <- signature_to_group[row_sigs]
    original_names <- names(aln)

  } else {
    # Utilize provided matrix and default mapping
    unique_bin <- topomatrix
    assigned_groups <- rownames(unique_bin)
    original_names <- rownames(unique_bin)
  }

  # Compute the pairwise distance matrix based on the selected evolutionary metric
  if (metric == "gaps") {
    dist_mat <- stats::dist(unique_bin, method = "binary")
  } else if (metric == "residues") {
    dist_mat <- stats::dist(1 - unique_bin, method = "binary")
  } else if (metric == "both") {
    dist_mat <- stats::dist(unique_bin, method = "manhattan") / ncol(unique_bin)
  } else {
    stop("Metric must be 'gaps', 'residues', or 'both'.")
  }

  # Attach underlying matrix and group metadata for seamless pipeline integration
  attr(dist_mat, "matrix") <- unique_bin
  attr(dist_mat, "metric") <- metric
  attr(dist_mat, "assigned_groups") <- assigned_groups
  attr(dist_mat, "original_names") <- original_names

  return(dist_mat)
}
