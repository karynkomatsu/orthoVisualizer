# Check that dataframe is formatted correctly
test_that("correct number of rows (ie. number of sequences)", {
  DNA_seq_path <-  system.file("extdata", "DNA_seq.fasta", package = "orthoVisualizer")
  expect_equal(nrow(quantSeq(DNA_seq_path, "CGCG")), 7)
})
test_that("correct number of columns", {
  DNA_seq_path <-  system.file("extdata", "DNA_seq.fasta", package = "orthoVisualizer")
  expect_equal(ncol(quantSeq(DNA_seq_path, "CGCG")), 5)
})

# Check that the calculations were done properly
test_that("occurrence quantification works for first row; CGCG", {
  DNA_seq_path <-  system.file("extdata", "DNA_seq.fasta", package = "orthoVisualizer")
  expect_equal((quantSeq(DNA_seq_path, "CGCG"))$numMotif[1], 3)
})
test_that("ratio caculation works for first row; CGCG", {
  DNA_seq_path <-  system.file("extdata", "DNA_seq.fasta", package = "orthoVisualizer")
  expect_equal((quantSeq(DNA_seq_path, "CGCG"))$ratioMotif[1], 0.02264151)
})
test_that("occurrence quantification works for first row; ATCGCG", {
  DNA_seq_path <-  system.file("extdata", "DNA_seq.fasta", package = "orthoVisualizer")
  expect_equal((quantSeq(DNA_seq_path, "ATCGCG"))$numMotif[1], 0)
})
test_that("ratio caculation works for first row; ATCGCG", {
  DNA_seq_path <-  system.file("extdata", "DNA_seq.fasta", package = "orthoVisualizer")
  expect_equal((quantSeq(DNA_seq_path, "ATCGCG"))$ratioMotif[1], 0)
})
test_that("occurrence quantification works for last row; CGCG", {
  DNA_seq_path <-  system.file("extdata", "DNA_seq.fasta", package = "orthoVisualizer")
  expect_equal((quantSeq(DNA_seq_path, "CGCG"))$numMotif[7], 19)
})
test_that("ratio caculation works for last row; CGCG", {
  DNA_seq_path <-  system.file("extdata", "DNA_seq.fasta", package = "orthoVisualizer")
  expect_equal((quantSeq(DNA_seq_path, "CGCG"))$ratioMotif[7], 0.05914397)
})
test_that("occurrence quantification works for last row; ATCGCG", {
  DNA_seq_path <-  system.file("extdata", "DNA_seq.fasta", package = "orthoVisualizer")
  expect_equal((quantSeq(DNA_seq_path, "ATCGCG"))$numMotif[7], 2)
})
test_that("ratio caculation works for last row; ATCGCG", {
  DNA_seq_path <-  system.file("extdata", "DNA_seq.fasta", package = "orthoVisualizer")
  expect_equal((quantSeq(DNA_seq_path, "ATCGCG"))$ratioMotif[7], 0.009338521)
})
