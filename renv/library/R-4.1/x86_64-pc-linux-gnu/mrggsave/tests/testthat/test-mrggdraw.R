library(testthat)
library(mrggsave)
library(ggplot2)

testthat::context("test-mrggsave-list")

set.seed(1100022)
data <- data.frame(x = rnorm(100), y = rnorm(100))
Script <- "test-mrggsave"
options(mrggsave.dir = tempdir())

p <- ggplot(data, aes(x = x, y = y)) + geom_point()

test_that("save a list", {
  l <- list(p,p,p)
  x <- mrggdraw(l)
  expect_length(x,3)
})

