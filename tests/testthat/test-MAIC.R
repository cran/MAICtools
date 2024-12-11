
test_that("Test MAIC functions", {

  require(tidyverse, quietly = TRUE, warn.conflicts = FALSE)

  data(IPD, package = "MAICtools")
  data(AgD_bl, package = "MAICtools")

  ## Baseline match-------------------------------------------------------------
  cov <- list(
    c("ECOG", "SMK", "METBRAIN"),
    c("BMI", "DIAG")
  )

  cov_all <- list(
    c("SEX", "ECOG", "SMK", "METBRAIN", "METLIVER"),
    c("BMI", "DIAG", "WEIGHT", "HEIGHT")
  )

  ##############################################################################
  pts <- estimate_weights(
    ipds = IPD, agds = AgD_bl, matching.list = cov,
    intervention.arm = TRT,
    comparator = STUDY, comparator.study = "Study XX-1", comparator.arm = TRT)

  missing_wt <- sum(is.na(pts$wt))

  expect_equal(missing_wt, 0)

  ##############################################################################
  ess <- estimate_ess(
    ipds_wts = pts, agds = AgD_bl,
    intervention.arm = TRT,
    comparator = STUDY, comparator.study = "Study XX-1", comparator.arm = TRT,
    comparator.n = N)

  ess_total <- ess[ess$TRT == "total", ]
  expect_true(all(ess_total$intervention_ess == 224))

  ##############################################################################
  weight_summary <- summarize_weights(ipds_wts = pts, intervention.arm = TRT)

  weight_active <- weight_summary[weight_summary$TRT == "active"
                                  & weight_summary$var == "Weights", ]
  weight_control <- weight_summary[weight_summary$TRT == "control"
                                   & weight_summary$var == "Weights", ]

  expect_equal(as.numeric(weight_active$mean), 0.66, tolerance = 0.01)
  expect_equal(as.numeric(weight_control$mean), 0.50, tolerance = 0.01)

})
