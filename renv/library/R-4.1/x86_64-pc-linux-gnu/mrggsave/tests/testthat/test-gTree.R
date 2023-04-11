library(testthat)
library(mrggsave)

testthat::context("gTree objects")

test_that("gTree", {
  pg <- grid::gTree(ggplot2::ggplot())
  expect_is(pg,"gTree")
  out <- mrggsave(pg, "foo.R", stem = "gTree", dir = tempdir())
  expect_equal(basename(out), "gTree.pdf")
})

