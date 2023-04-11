library(testthat)
library(mrggsave)
library(ggplot2)
library(grid)

testthat::context("test-devices")

set.seed(1100022)
data <- data.frame(x = rnorm(100), y = rnorm(100))
Script <- "test-mrggsave"
options(mrggsave.dir = tempdir(), mrg.script = "test.R")

pg <- ggplot(data, aes(x = x, y = y)) + geom_point()


test_that("png", {
  foo <- mrggsave(pg, stem = "foo", dev="png")
  expect_equal(basename(foo), "foo.png")
})

test_that("bmp", {
  foo <- mrggsave(pg, stem = "foo", dev="bmp")
  expect_equal(basename(foo), "foo.bmp")
})

test_that("eps", {
  foo <- mrggsave(pg, stem = "foo", dev="eps")
  expect_equal(basename(foo), "foo.eps")
})

test_that("ps", {
  foo <- mrggsave(pg, stem = "foo", dev="ps")
  expect_equal(basename(foo), "foo.ps")
})

test_that("tiff", {
  foo <- mrggsave(pg, stem = "foo", dev="tiff")
  expect_equal(basename(foo), "foo.tiff")
})

test_that("pdf", {
  foo <- mrggsave(pg, stem = "foo", dev="pdf")
  expect_equal(basename(foo), "foo.pdf")
})

test_that("jpeg", {
  foo <- mrggsave(pg, stem = "foo", dev="jpeg")
  expect_equal(basename(foo), "foo.jpeg")
})

test_that("svg", {
  foo <- mrggsave(pg, stem = "foo", dev="svg")
  expect_equal(basename(foo), "foo.svg")
})

test_that("cairo_pdf", {
  foo <- mrggsave(pg, stem = "foo", dev="cairo_pdf")
  expect_equal(basename(foo), "foo.pdf")
})

test_that("save multiple plots to single file with cairo_pdf", {
  foo <- mrggsave(list(pg,pg,pg), stem = "multi-cairo", dev = "cairo_pdf")
  expect_identical(basename(foo), "multi-cairo.pdf")
})

test_that("save to multiple devices", {
  foo <- mrggsave(list(pg,pg), stem = "multiple", dev = "pdf,png")
  expect_true(file.exists(foo[[1]]))
  expect_true(file.exists(foo[[2]]))
  expect_true(file.exists(foo[[3]]))
  foo <- basename(foo)
  expect_equal(foo, c("multiple.pdf", "multiple001.png", "multiple002.png"))
})
