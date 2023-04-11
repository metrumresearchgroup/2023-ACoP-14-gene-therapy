library(testthat)
library(mrggsave)

testthat::context("test-ggsurv")

test_that("ggsurvplot", {
  skip_if_not(requireNamespace("survminer",quietly=TRUE))

  data(cancer, package = "survival")

  fit <- survival::survfit(survival::Surv(time, status) ~ sex, data = cancer)

  # With a new theme
  pp1 <- survminer::ggsurvplot(fit, data = cancer, ggtheme = ggplot2::theme_bw())
  x <- mrggsave(pp1, script = "test-ggsurv.R", dir = tempdir(), tag = "with_theme")
  expect_identical(basename(x),"test-ggsurv-with_theme.pdf")

  # With the default theme
  pp2 <- survminer::ggsurvplot(fit, data = cancer)
  x <- mrggsave(pp2, script = "test-ggsurv.R", dir = tempdir(), tag = "no_theme")
  expect_identical(basename(x),"test-ggsurv-no_theme.pdf")
})
