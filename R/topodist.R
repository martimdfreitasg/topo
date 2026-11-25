#' Calculate Alignment Topological Distances
#'
#' Computes pairwise topological distances between sequences in a multiple
#' sequence alignment based on binary gap/residue patterns.
#'
#' @details
#' Each aligned sequence is encoded as a binary vector where gaps (`"-"`) are
#' represented as `1` and residue positions as `0`. Sequences sharing
#' the same pattern are collapsed into a single
#' representative topology before distance computation to avoid redundant
#' calculations and zero-length clusters; the mapping between original headers
#' and representative groups is stored in attributes of the returned object.
#'
#' Three complementary distance metrics are supported:
#' \itemize{
#'   \item `"gaps"`: Asymmetric Jaccard distance computed on gap positions
#'     (`stats::dist(x, method = "binary")`). Shared residue positions
#'     (double-zeros) are ignored, focusing on shared insertion/deletion
#'     (indels).
#'   \item `"residues"`: Asymmetric Jaccard distance computed on
#'     residue positions (`stats::dist(1 - x, method = "binary")`). Shared gap
#'     blocks are ignored, which is suited to sparse alignments containing
#'     lineage-specific insertions or partial sequences.
#'   \item `"both"`: Normalized Hamming distance
#'     (`stats::dist(x, method = "manhattan") / ncol(x)`), weighting shared gaps
#'     and residues across all alignment columns.
#' }
#'
#' @param fasta_path Path to a FASTA alignment file, or an `ape::AAbin` /
#'   `ape::DNAbin` object. Optional if `topomatrix` is provided.
#' @param topomatrix Optional numeric binary matrix (`0` = residue, `1` = gap)
#'   used directly instead of reading `fasta_path` (primarily for internal
#'   bootstrapping or pre-processed matrices).
#' @param metric Character string specifying the topological distance metric:
#'   one of `"gaps"` (default), `"residues"`, or `"both"`.
#'
#' @return A `stats::dist` object of pairwise topological distances among unique
#'   sequence topologies, with the following attributes attached:
#'   \itemize{
#'     \item `matrix`: Deduplicated binary matrix (`0` = residue, `1` = gap).
#'     \item `metric`: Distance metric used.
#'     \item `assigned_groups`: Character vector mapping each input sequence to
#'       its representative tip/group label.
#'     \item `original_names`: Character vector of original sequence headers.
#'   }
#'
#' @references
#' Paradis, E., & Schliep, K. (2019). ape 5.0: an environment for modern
#' phylogenetics and evolutionary analyses in R. *Bioinformatics*, 35(3),
#' 526--528. \doi{10.1093/bioinformatics/bty633}
#'
#' @export
#' @examples
#' fasta_path = system.file("extdata", "test.fasta", package = "topo")
#' gap_dist = topodist(fasta_path, metric = "gaps")
#' gap_dist
topodist = function(fasta_path = NULL,
                     topomatrix = NULL,
                     metric = c("gaps", "residues", "both")) {

  metric = match.arg(metric)

  if (is.null(topomatrix)) {
    if (is.null(fasta_path)) {
      stop("Either 'fasta_path' or 'topomatrix' must be provided.", call. = FALSE)
    }

    # Parse supported alignment inputs into a named list of single-character vectors
    if (is.character(fasta_path)) {
      if (length(fasta_path) != 1L || !nzchar(fasta_path) || !file.exists(fasta_path)) {
        stop("FASTA file not found: ", paste(fasta_path, collapse = ", "), call. = FALSE)
      }
      aln_raw = ape::read.FASTA(fasta_path, type = "AA")
      if (length(aln_raw) == 0L) {
        stop("No sequences found in the FASTA file.", call. = FALSE)
      }
      aln = as.character(as.list(aln_raw))
    } else if (inherits(fasta_path, c("AAbin", "DNAbin"))) {
      aln = as.character(as.list(fasta_path))
      if (length(aln) == 0L) {
        stop("Input alignment object contains no sequences.", call. = FALSE)
      }
    } else {
      stop(
        "'fasta_path' must be a path to a FASTA file or an ape 'AAbin'/'DNAbin' object.",
        call. = FALSE
      )
    }

    if (length(aln) < 2L) {
      stop("At least 2 sequences are required to compute topological distances.", call. = FALSE)
    }

    original_names = names(aln)
    if (is.null(original_names)) {
      original_names = paste0("Seq_", seq_along(aln))
      names(aln) = original_names
    }
    if (anyDuplicated(original_names)) {
      dups = unique(original_names[duplicated(original_names)])
      stop(
        "Duplicate sequence headers detected: ",
        paste(utils::head(dups, 5L), collapse = ", "),
        if (length(dups) > 5L) ", ..." else "",
        ". All sequence headers must be unique.",
        call. = FALSE
      )
    }

    # Pad uneven sequence lengths with terminal gaps for matrix conversion
    seq_lens = lengths(aln)
    max_len = max(seq_lens)
    if (max_len == 0L) {
      stop("Input sequences have zero length.", call. = FALSE)
    }
    if (length(unique(seq_lens)) > 1L) {
      aln = lapply(aln, function(x) c(x, rep("-", max_len - length(x))))
    }

    # Coerce alignment to a binary matrix (1 = Gap, 0 = Residue)
    bin_mat = (do.call(rbind, aln) == "-") + 0
    rownames(bin_mat) = original_names

    # Deduplicate identical topologies into representative groups
    row_sigs = apply(bin_mat, 1L, paste0, collapse = "")
    unique_idx = !duplicated(row_sigs)
    unique_bin = bin_mat[unique_idx, , drop = FALSE]

    if (nrow(unique_bin) < 2L) {
      stop(
        "All sequences share the exact same gap/residue topology (only 1 unique pattern found). ",
        "At least 2 distinct topologies are required for distance calculation and clustering.",
        call. = FALSE
      )
    }

    unique_sigs = row_sigs[unique_idx]
    sig_id = match(row_sigs, unique_sigs)
    sig_counts = tabulate(sig_id, nbins = length(unique_sigs))
    is_group = sig_counts > 1L

    # Assign collision-free "Group k" labels to duplicated topologies
    if (any(is_group)) {
      existing_singletons = rownames(unique_bin)[!is_group]
      n_groups = sum(is_group)
      group_names = character(n_groups)
      k = 1L
      for (g in seq_len(n_groups)) {
        candidate = paste0("Group ", k)
        while (candidate %in% existing_singletons) {
          k = k + 1L
          candidate = paste0("Group ", k)
        }
        group_names[g] = candidate
        k = k + 1L
      }
      rownames(unique_bin)[is_group] = group_names
    }

    assigned_groups = unname(rownames(unique_bin)[sig_id])

  } else {
    if (!is.matrix(topomatrix) || !is.numeric(topomatrix)) {
      stop("'topomatrix' must be a numeric binary matrix (0 = residue, 1 = gap).", call. = FALSE)
    }
    if (anyNA(topomatrix) || !all(topomatrix %in% c(0, 1))) {
      stop("'topomatrix' must contain only binary values (0 and 1) with no missing values.", call. = FALSE)
    }
    if (nrow(topomatrix) < 2L) {
      stop("'topomatrix' must contain at least 2 rows.", call. = FALSE)
    }
    if (ncol(topomatrix) < 1L) {
      stop("'topomatrix' must contain at least 1 column.", call. = FALSE)
    }

    unique_bin = topomatrix
    if (is.null(rownames(unique_bin))) {
      rownames(unique_bin) = paste0("Seq_", seq_len(nrow(unique_bin)))
    }
    assigned_groups = rownames(unique_bin)
    original_names = rownames(unique_bin)
  }

  # Compute pairwise distances
  dist_mat = switch(
    metric,
    gaps = stats::dist(unique_bin, method = "binary"),
    residues = stats::dist(1 - unique_bin, method = "binary"),
    both = stats::dist(unique_bin, method = "manhattan") / ncol(unique_bin)
  )

  attr(dist_mat, "matrix") = unique_bin
  attr(dist_mat, "metric") = metric
  attr(dist_mat, "assigned_groups") = assigned_groups
  attr(dist_mat, "original_names") = original_names

  return(dist_mat)
}
