library(testthat)
library(mrggsave)
library(ggplot2)

testthat::context("test-mrggpage")

set.seed(1100022)
data <- data.frame(x = rnorm(10), y = rnorm(10))
p <- ggplot(data = data, aes(x,y)) + geom_point()
Script <- "test-mrggsave"
options(mrggsave.dir = tempdir())

test_that("plots are arranged on a single page", {
  x <- list(p,p,p,p)
  y <- mrggpage(x, ncol = 2)
  expect_is(y, "gtable")
  expect_equal(dim(y), c(2,2))

})

test_that("plots are arranged on a multiple pages", {
  x <- list(p,p,p,p)
  y <- mrggpage(x, ncol = 2, nrow = 1, multiple = TRUE)
  expect_is(y, "arrangelist")
  expect_equal(length(y), 2)
  expect_is(y[[1]], "gtable")
  expect_equal(dim(y[[1]]), c(1,2))
})

