test_that("topodist computes valid distances across all metrics and validates metric early", {
  fasta_path <- system.file("extdata", "test.fasta", package = "topo")

  for (m in c("gaps", "residues", "both")) {
    d <- topodist(fasta_path, metric = m)
    expect_s3_class(d, "dist")
    expect_equal(attr(d, "metric"), m)
    expect_true(is.matrix(attr(d, "matrix")))
    expect_false(anyNA(d))
  }

  expect_error(topodist("nonexistent_file.fasta", metric = "invalid"), "should be one of")
  expect_error(topodist("nonexistent_file.fasta", metric = "gaps"), "FASTA file not found")
})

test_that("topodist handles edge cases and ape AAbin objects cleanly", {
  # Identical topologies
  f_same <- tempfile(fileext = ".fasta")
  writeLines(c(">seq1", "AC-GT", ">seq2", "TG-CA"), f_same)
  on.exit(unlink(f_same), add = TRUE)
  expect_error(topodist(f_same), "exact same gap/residue topology")

  # Duplicate headers
  f_dup <- tempfile(fileext = ".fasta")
  writeLines(c(">seq1", "AC-GT", ">seq1", "A--GT"), f_dup)
  on.exit(unlink(f_dup), add = TRUE)
  expect_error(topodist(f_dup), "Duplicate sequence headers")

  # Collision with an existing "Group 1" singleton header
  f_grp <- tempfile(fileext = ".fasta")
  writeLines(c(">Group 1", "A---T", ">seqA", "AC-GT", ">seqB", "TG-CA"), f_grp)
  on.exit(unlink(f_grp), add = TRUE)
  d_grp <- topodist(f_grp, metric = "gaps")
  expect_equal(anyDuplicated(rownames(attr(d_grp, "matrix"))), 0L)
  expect_true("Group 1" %in% rownames(attr(d_grp, "matrix")))
  expect_true("Group 2" %in% rownames(attr(d_grp, "matrix")))

  # Direct ape::AAbin input
  fasta_path <- system.file("extdata", "test.fasta", package = "topo")
  aa_obj <- ape::read.FASTA(fasta_path, type = "AA")
  d_aa <- topodist(aa_obj, metric = "gaps")
  expect_s3_class(d_aa, "dist")
})

test_that("topoclust exposes deepSplit and minClusterSize", {
  fasta_path <- system.file("extdata", "test.fasta", package = "topo")
  set.seed(123)
  cl <- topoclust(fasta_path, metric = "gaps", deepSplit = 4, minClusterSize = 2)
  expect_s3_class(cl, "data.frame")
  expect_named(cl, c("Original_Sequence", "Assigned_Group", "Cluster"))
  expect_error(topoclust(fasta_path, deepSplit = 9), "deepSplit")
  expect_error(topoclust(fasta_path, minClusterSize = 0), "minClusterSize")
})

test_that("topodend uses percentage bootstrap, respects metric_used, and preserves ultrametricity", {
  fasta_path <- system.file("extdata", "test.fasta", package = "topo")
  d_res <- topodist(fasta_path, metric = "residues")

  set.seed(42)
  dend <- topodend(dist_mat = d_res, bootstrap = 20, boot_threshold = 70)
  expect_s3_class(dend$as_phylo, "phylo")
  expect_s3_class(dend$as_ggtree, "ggplot")

  # All tip depths remain ultrametric after collapsing branches
  depths <- ape::node.depth.edgelength(dend$as_phylo)[seq_along(dend$as_phylo$tip.label)]
  expect_equal(max(depths) - min(depths), 0, tolerance = 1e-8)

  # Numeric node labels must be on the 0-100 percentage scale
  non_empty <- dend$as_phylo$node.label[nzchar(dend$as_phylo$node.label)]
  if (length(non_empty) > 0L) {
    vals <- as.numeric(non_empty)
    expect_true(all(vals >= 70 & vals <= 100))
  }
})

test_that("topomap and topoplot validate highlights and render without ggplot2 deprecation warnings", {
  fasta_path <- system.file("extdata", "test.fasta", package = "topo")

  expect_error(
    topomap(fasta_path, highlight_regions = list(1:999999)),
    "out-of-bounds"
  )

  expect_warning(
    p_ov <- topomap(fasta_path, highlight_regions = list(1:10, 5:15)),
    "overlaps"
  )
  expect_s3_class(p_ov, "ggplot")

  expect_no_warning(
    p_full <- topoplot(fasta_path, metric = "gaps")
  )
  expect_s3_class(p_full, "ggplot")
})

test_that("topodist rejects NA and non-binary values in topomatrix", {
  bad_na <- matrix(c(0, 1, NA, 0), nrow = 2)
  bad_val <- matrix(c(0, 2, 1, 0), nrow = 2)

  expect_error(topodist(topomatrix = bad_na), "binary values")
  expect_error(topodist(topomatrix = bad_val), "binary values")
})

test_that("topodend anchors distance labels at 0.00 in full pairwise scale without messages", {
  m <- matrix(c(1, 0, 0, 1), nrow = 2, dimnames = list(c("A", "B"), NULL))
  d <- topodist(topomatrix = m, metric = "gaps") # Jaccard distance = 1.0

  expect_no_message(dend <- topodend(dist_mat = d))
  built <- ggplot2::ggplot_build(dend$as_ggtree)
  all_labels <- unlist(lapply(built$data, function(df) {
    if ("label" %in% names(df)) as.character(df$label)
  }))

  expect_true("0.00" %in% all_labels)
  expect_true("1.00" %in% all_labels)
})
