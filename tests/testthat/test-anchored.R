
test_that("Test anchored functions", {

  require(tidyverse, quietly = TRUE, warn.conflicts = FALSE)

  data(pts, package = "MAICtools")
  data(AgD_eff, package = "MAICtools")

  pts <- pts %>%
    mutate(TRT = factor(TRT, levels = c("control", "active")))

  ##############################################################################
  results1 <- anchored_maic(
    ipds_wts = pts, intervention.arm = TRT,
    agds_eff = AgD_eff, comparator = STUDY,
    comparator.study = "Study XX-1",
    ipds.param.var = PARAMCD, ipds.param = "OS",
    agds.param.var = PARAM, agds.param = "OS",
    agds.estimate = EST, agds.ci.lower = CIL, agds.ci.upper = CIU,
    time = AVAL, status = CNSR, event = 0,
    stralist = "BPDL1, CNSBRAIN, AGEGR", dtype = "HR",
    wt.col = wt, CIw = 0.95, digits = 2)

  results1_unadjusted <- results1[results1$outcome == "unadjusted", ]
  results1_weighted <- results1[results1$outcome == "weighted", ]

  expect_equal(as.numeric(results1_unadjusted$effects), 0.81, tolerance = 0.01)
  expect_equal(as.numeric(results1_weighted$effects), 0.78, tolerance = 0.01)

  ##############################################################################

  results2 <- anchored_maic(
    ipds_wts = pts, intervention.arm = TRT,
    agds_eff = AgD_eff, comparator = STUDY,
    comparator.study = "Study XX-1",
    agds.param.var = PARAM, agds.param = "ORR",
    agds.estimate = EST, agds.ci.lower = CIL, agds.ci.upper = CIU,
    response = RESP,
    stralist = "BPDL1, CNSBRAIN, AGEGR", dtype = "OR",
    wt.col = wt, CIw = 0.95, digits = 2)

  results2_unadjusted <- results2[results2$outcome == "unadjusted", ]
  results2_weighted <- results2[results2$outcome == "weighted", ]

  expect_equal(as.numeric(results2_unadjusted$effects), 0.80, tolerance = 0.01)
  expect_equal(as.numeric(results2_weighted$effects), 0.82, tolerance = 0.01)

})
