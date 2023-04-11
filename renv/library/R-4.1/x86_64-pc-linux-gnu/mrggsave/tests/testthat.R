Sys.setenv("R_TESTS" = "")
library(testthat)
library(mrggsave)

test_check("mrggsave", reporter="summary")

