
test_that("Test anchored functions", {

  require(tidyverse, quietly = TRUE, warn.conflicts = FALSE)

  data(unpts, package = "MAICtools")

  ##############################################################################
  results3 <- unanchored_maic(
    unds_wts = unpts, unds.arm = ARM,
    comparator.study = "IMpower133",
    unds.param.var = PARAMCD, unds.param = "OS",
    time = AVAL, status = CNSR, event = 0,
    dtype = "HR")

  results3_unadjusted <- results3[results3$outcome == "unadjusted", ]
  results3_weighted <- results3[results3$outcome == "weighted", ]

  expect_equal(as.numeric(results3_unadjusted$effects), 0.75, tolerance = 0.01)
  expect_equal(as.numeric(results3_weighted$effects), 0.66, tolerance = 0.01)

  ##############################################################################

  results4 <- unanchored_maic(
    unds_wts = unpts, unds.arm = ARM,
    # unds.param = "ORR",
    # comparator.study = "IMpower133",
    response = CNSR,
    dtype = "OR")

  results4_unadjusted <- results4[results4$outcome == "unadjusted", ]
  results4_weighted <- results4[results4$outcome == "weighted", ]

  expect_equal(as.numeric(results4_unadjusted$effects), 1.84, tolerance = 0.01)
  expect_equal(as.numeric(results4_weighted$effects), 2.14, tolerance = 0.01)

})
