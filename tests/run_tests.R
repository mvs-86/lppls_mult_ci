setwd("C:/Users/marcus/Projetos/lppls_mult_ci")
library(testthat)

cat("--- test_lppls_fit.R ---\n")
test_file("tests/test_lppls_fit.R")

cat("\n--- test_lppls_ci.R ---\n")
test_file("tests/test_lppls_ci.R")
