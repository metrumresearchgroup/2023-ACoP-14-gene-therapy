library(testthat)
library(mrggsave)
library(ggplot2)
library(glue)

testthat::context("test-filename")

data <- data.frame(x = c(1,2,3), y = c(4,5,6))

p  <- ggplot(data) + geom_point(aes(x,y))
p1 <- ggplot(data) + geom_point(aes(x,y))
p2 <- ggplot(data) + geom_point(aes(x,y))
p3 <- ggplot(data) + geom_point(aes(x,y))

options(mrggsave.dir = normalizePath(tempdir()))

test_that("variable gets glued into stem", {
  runn <- 1234
  ans <- mrggsave(p, script = "test-filename.R", stem = "save_{runn}", dir = tempdir())
  expect_equal(basename(ans), "save_1234.pdf")
})

test_that("variable gets glued into tag", {
  runn <- 1234
  ans <- mrggsave(p, script = "file.R", tag = "save_{runn}")
  expect_equal(basename(ans), "file-save_1234.pdf")
})

test_that("an environment gets passed to glue", {
  env <- list(runn=4321)
  ans <- mrggsave(p, script = "file.R", tag = "save_{runn}", envir = env)
  expect_equal(basename(ans), "file-save_4321.pdf")
})

test_that("vector stem gets collapsed", {
  ans <- mrggsave(p, script = "test-filename.R", stem = c("a", 101, "b"))
  expect_equal(basename(ans), "a-101-b.pdf")
})

test_that("vector tag gets collapsed", {
  ans <- mrggsave(p, script = "test-filename.R", tag = c("a", 101, "b"))
  expect_equal(basename(ans), "test-filename-a-101-b.pdf")
})

test_that("plots get named by object", {
  p1 <- p2 <- p3 <- p
  l <- named_plots(p1,p2,p3, tag = "bbb")
  expect_identical(names(l), c("p1-bbb", "p2-bbb", "p3-bbb"))
  expect_length(l,3)
  cl <- sapply(l, is.ggplot)
  expect_true(all(cl))
})

test_that("change file name separator", {
  out1 <- basename(mrggsave(p1, tag = "1", script = "foo.R", dev = "bmp"))
  expect_equal(out1, "foo-1.bmp")
  mrggsave:::output_file_sep("_")
  out2 <- basename(mrggsave(p1, tag = "1", script = "foo.R", dev = "bmp"))
  expect_equal(out2, "foo_1.bmp")
  mrggsave:::output_file_sep()
})

test_that("named_plots returns an object with class", {
  p1 <- p2 <- p3 <- p
  ans <- named_plots(p1,p2,p3)
  expect_is(ans, "named-plots")
  ans <- named_plots(p1, add_context = TRUE)
  expect_is(ans, "named-plots")
  expect_is(ans, "needs-context")
})

test_that("named_plots input auto uses names", {
  p1 <- p2 <- p3 <- p
  inpt <- named_plots(a = p1, p2, ggplot(mtcars))
  expect_equal(names(inpt), c("a", "p2", "ggplot"))
  ans <- mrggsave(inpt, script = "test-filename.R")
  ans <- basename(ans)
  expect_equal(ans, c("a.pdf", "p2.pdf", "ggplot.pdf"))
  inpt <- named_plots(p1, add_context = TRUE)
  ans <- mrggsave(inpt, script = "scrname.R")
  ans <- basename(ans)
  expect_equal(ans, "scrname-p1.pdf")
})

test_that("named_plots names are sanitized", {
  dv_pred <- p1
  inpt <- named_plots(dv_pred, "a b.-c" = p2)
  expect_equal(names(inpt), c("dv-pred", "a-b-c"))
})

options(mrggsave.file.tolower = TRUE)
test_that("option to make lower", {
  inpt <- named_plots(EV.PREd = p1)
  ans <- mrggsave(inpt, script = "test-filename.R")
  expect_equal(basename(ans), "ev-pred.pdf")
  ans <- mrggsave(p1, script = "test-filename", stem = "ABCDE")
  expect_equal(basename(ans), "abcde.pdf")
})
options(mrggsave.file.tolower = NULL)

test_that("passing a named list with use_names set to TRUE", {
  p1 <- p2 <- p3 <- p
  ans <- mrggsave(list(a = p1, b = p1), script = "blah.R", use_names = TRUE)
  expect_equal(basename(ans), c("a.pdf","b.pdf"))
  ans <- mrggsave(list(a = p1, b = p1), script = "blah.R", use_names = FALSE)
  def_stem <- formals(mrggsave_common)$stem
  check <- paste0(def_stem, ".pdf")
  expect_equal(basename(ans), check)
})

test_that("glue file name within a function", {
  fun <- function(x, runno) {
    mrggsave(p, stem = "plot-{runno}", dir = tempdir(), script = "test.R")
  }
  p1 <- p
  ans <- fun(p1, "112")
  expect_equal(basename(ans), "plot-112.pdf")
})

test_that("glue file name with mrggsave_last", {
  p <- ggplot(data) + geom_point(aes(x,y))
  runno <- 5678
  ans <- mrggsave_last(
    stem = "plot-{runno}",
    dir = tempdir(),
    script = "test.R"
  )
  expect_equal(basename(ans), "plot-5678.pdf")
})
