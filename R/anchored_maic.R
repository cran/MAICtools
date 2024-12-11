#' @title Conduct Anchored Matching-Adjusted Indirect Comparison (MAIC).
#'
#' @description
#' The endpoint of interest is either time-to-event (e.g., overall survival)
#' or binary (e.g., objective tumor response). The methods described in this
#' documentation are based on those originally outlined by Signorovitch et al.,
#' 2012, and further detailed in the National Institute for Health and Care
#' Excellence (NICE) Decision Support Unit (DSU) Technical Support Document
#' (TSD) 18.
#'
#' @param ipds_wts A data frame containing individual patient data from the
#' intervention study, with a column containing the estimated weights (derived
#' using \code{\link{estimate_weights}}).
#' @param intervention.arm The name of the grouping column in the data frame
#' specified by `ipds_wts`, e.g., intervention.arm = TRT. The default is `TRT`.
#' @param agds_eff A data frame containing aggregate efficacy results from the
#' comparator study.
#' @param comparator The name of the study column in the data frame specified
#' by `agds_eff`, e.g., `comparator = STUDY`. The default is `STUDY`.
#' @param comparator.study A character specifying the comparator study, which
#' must be quoted and cannot be empty (e.g., `comparator.study = "Study XX-1"`).
#' This is the value of the study column in `agds_eff` set by the `comparator`
#' parameter.
#' @param ipds.param.var The name of the column that specifies only a subset
#' of the `ipds_wts` to be used.
#' @param ipds.param A character specifying the subset of the rows to be used.
#' This is the value of the column set by the `ipds.param.var`.
#' @param agds.param.var The name of the column that specifies only a specific
#' result of the `agds_eff` to be used.
#' @param agds.param A character specifying the subset of the rows to be used.
#' This is the value of the column set by the `agds.param.var`.
#' @param agds.estimate The column name of the point estimate of the effect
#' size.
#' @param agds.ci.lower The column name for the lower confidence limit of the
#' point estimate of the effect size.
#' @param agds.ci.upper The column name for the upper confidence limit of the
#' point estimate of the effect size.
#' @param time The name of the survival or follow-up time column in the
#' `ipds_wts`.
#' @param status The status indicator, normally 0 = event, 1 = censored. Can be
#' reset using the `event` parameter.
#' @param event A numeric value that represents the survival status, 0 = event,
#' 1 = censored.
#' @param response The name of the response status column in the `ipds_wts`.
#' @param stralist A string specifying the stratification factors in a
#' stratified analysis, e.g., `stralist = "BPDL1, CNSBRAIN, AGEGR"`.
#' @param dtype Two options are available: "HR" and "OR". The default is "HR".
#' @param wt.col The name of the estimated weights column in the data frame
#' specified by `ipds_wts`. The default is `wt`.
#' @param CIw The numeric value specifying the width of the confidence
#' interval, with a default of 0.95.
#' @param digits Specify the number of decimal places for the output results.
#'
#' @importFrom rlang enquo as_name quo_is_null eval_tidy get_expr as_string quo
#' @importFrom assertthat assert_that
#' @importFrom dplyr transmute mutate select filter bind_rows
#' @importFrom stats qnorm coef glm binomial as.formula
#' @importFrom survival Surv strata coxph
#' @importFrom broom tidy
#' @importFrom stringr str_detect
#' @importFrom data.table :=
#'
#' @return A data frame containing the anchored matching-adjusted indirect
#' comparison results.
#' @export
#'
#' @examples
#' \donttest{
#' results1 <- anchored_maic(
#'   ipds_wts = pts, intervention.arm = TRT,
#'   agds_eff = AgD_eff, comparator = STUDY,
#'   comparator.study = "Study XX-1",
#'   ipds.param.var = PARAMCD, ipds.param = "OS",
#'   agds.param.var = PARAM, agds.param = "OS",
#'   agds.estimate = EST, agds.ci.lower = CIL, agds.ci.upper = CIU,
#'   time = AVAL, status = CNSR, event = 0,
#'   stralist = "BPDL1, CNSBRAIN, AGEGR", dtype = "HR",
#'   wt.col = wt, CIw = 0.95, digits = 2)
#'
#' results1
#'
#' results2 <- anchored_maic(
#'   ipds_wts = pts, intervention.arm = TRT,
#'   agds_eff = AgD_eff, comparator = STUDY,
#'   comparator.study = "Study XX-1",
#'   agds.param.var = PARAM, agds.param = "ORR",
#'   agds.estimate = EST, agds.ci.lower = CIL, agds.ci.upper = CIU,
#'   response = RESP,
#'   stralist = "BPDL1, CNSBRAIN, AGEGR", dtype = "OR",
#'   wt.col = wt, CIw = 0.95, digits = 2)
#'
#' results2
#' }
#'
#' @name anchored_maic
NULL

# 声明全局变量
utils::globalVariables(c("TRT", "STUDY", "wt", "log_coef", "log_coef_se",
                         "param", "type", "effects", "lower", "upper",
                         ".", "term", "estimate", "std.error"))

anchored_maic <- function(
    ipds_wts, intervention.arm = TRT,
    agds_eff, comparator = STUDY,
    comparator.study = NA,
    ipds.param.var = NULL, ipds.param = NULL,
    agds.param.var = NULL, agds.param = NULL,
    agds.estimate, agds.ci.lower, agds.ci.upper,
    time = NULL, status = NULL, event = 0, response = NULL,
    stralist = NULL, dtype = "HR",
    wt.col = wt, CIw = 0.95, digits = 2) {

  # 使用enquo捕获参数符号
  intervention.arm <- rlang::enquo(intervention.arm)
  comparator <- rlang::enquo(comparator)
  ipds.param.var <- rlang::enquo(ipds.param.var)
  agds.param.var <- rlang::enquo(agds.param.var)
  agds.estimate <- rlang::enquo(agds.estimate)
  agds.ci.lower <- rlang::enquo(agds.ci.lower)
  agds.ci.upper <- rlang::enquo(agds.ci.upper)
  time <- rlang::enquo(time)
  status <- rlang::enquo(status)
  response <- rlang::enquo(response)
  wt.col <- rlang::enquo(wt.col)

  # checking--------------------------------------------------------------------
  # 检查ipds和agds是否为data.frame类型
  assertthat::assert_that(is.data.frame(ipds_wts),
                          msg = "'ipds_wts' is expected to be a data frame")

  assertthat::assert_that(is.data.frame(agds_eff),
                          msg = "'agds_eff' is expected to be a data frame")


  # 检查变量intervention.arm、ipds.param.var、wt.col是否在数据集ipds_wts中
  assertthat::assert_that(rlang::as_name(intervention.arm) %in% names(ipds_wts),
                          msg = paste0(rlang::as_name(intervention.arm),
                                       " can not be found in ipds_wts"))
  assertthat::assert_that(rlang::as_name(wt.col) %in% names(ipds_wts),
                          msg = paste0(rlang::as_name(wt.col),
                                       " can not be found in ipds_wts"))

  # 当dtype = "HR"时，检查变量time、status是否在数据集ipds_wts中
  if (dtype == "HR") {
    assertthat::assert_that(rlang::as_name(time) %in% names(ipds_wts),
                            msg = paste0(rlang::as_name(time),
                                         " can not be found in ipds_wts"))
    assertthat::assert_that(rlang::as_name(status) %in% names(ipds_wts),
                            msg = paste0(rlang::as_name(status),
                                         " can not be found in ipds_wts"))
  }

  # 当dtype = "OR"时，检查变量response是否在数据集ipds_wts中
  if (dtype == "OR") {
    assertthat::assert_that(rlang::as_name(response) %in% names(ipds_wts),
                            msg = paste0(rlang::as_name(response),
                                         " can not be found in ipds_wts"))
  }

  # 检查变量comparator、agds.param.var是否在数据集agds_eff中
  assertthat::assert_that(rlang::as_name(comparator) %in% names(agds_eff),
                          msg = paste0(rlang::as_name(comparator),
                                       " can not be found in agds_eff"))
  assertthat::assert_that(rlang::as_name(agds.param.var) %in% names(agds_eff),
                          msg = paste0(rlang::as_name(agds.param.var),
                                       " can not be found in agds_eff"))

  # 必要的函数定义
  estimate_itc <- function(intervention = NULL, comparator = NULL, outcome) {

    # 默认值设置
    if (rlang::quo_is_null(ipds.param.var)) {
      ipds.param.var <- quo(NULL)
    }
    if (is.null(ipds.param)) {
      ipds.param <- " "  # 适当的默认参数值
    }


    lower_ci <- paste0("lower_", CIw, "ci")
    upper_ci <- paste0("upper_", CIw, "ci")

    itc <- data.frame(
      log_coef = intervention$log_coef - comparator$log_coef,
      log_coef_se = sqrt(intervention$log_coef_se^2 + comparator$log_coef_se^2)
    ) %>%
      dplyr::transmute(
        coef = exp(log_coef),
        lower = exp(log_coef - stats::qnorm(1 - (1 - CIw) / 2) * log_coef_se),
        upper = exp(log_coef + stats::qnorm(1 - (1 - CIw) / 2) * log_coef_se)
      ) %>%
      dplyr::mutate(
        comparator = comparator.study,
        param = ipds.param,
        type = "anchored",
        outcome = outcome
      ) %>%
      dplyr::mutate(
        effects = sprintf(paste0("%.", digits, "f"), coef),
        !!lower_ci := sprintf(paste0("%.", digits, "f"), lower),
        !!upper_ci := sprintf(paste0("%.", digits, "f"), upper)
      ) %>%
      dplyr::select(
        comparator, param, type, outcome, effects, !!lower_ci, !!upper_ci
      )

    return(itc)
  }

  if (dtype == "HR") {

    # comparator study
    comparator_efficacy <- agds_eff %>%
      dplyr::filter(rlang::eval_tidy(comparator, .) == comparator.study &
                      rlang::eval_tidy(agds.param.var, .) == agds.param) %>%
      dplyr::transmute(
        log_coef = log(!!agds.estimate),
        log_coef_se = (log(!!agds.ci.upper)
                       - log(!!agds.ci.lower)) /
          (2 * stats::qnorm(1 - (1 - CIw) / 2))
      )

    # intervention study to obtain the HR (Hazard Ratio)
    formula <- stats::as.formula(
      paste0("survival::Surv(", rlang::as_string(rlang::get_expr(time)), ", ",
             rlang::as_string(rlang::get_expr(status)), "== ", event, ") ~",
             rlang::as_string(rlang::get_expr(intervention.arm)),
             if (!is.null(stralist)) paste0(" + strata(", stralist, ")")
             else ""))

    intervention_cox <- survival::coxph(
      formula, data = ipds_wts,
      subset = if (!is.null(ipds.param) && !is.null(ipds.param.var)) {
        rlang::eval_tidy(rlang::enquo(ipds.param.var), ipds_wts) == ipds.param
      } else {
        TRUE
      },
      method = "efron")

    intervention_cox_wtd <- survival::coxph(
      formula, data = ipds_wts,
      subset = if (!is.null(ipds.param) && !is.null(ipds.param.var)) {
        rlang::eval_tidy(rlang::enquo(ipds.param.var), ipds_wts) == ipds.param
      } else {
        TRUE
      },
      method = "efron",
      weights = ipds_wts[[rlang::as_name(wt.col)]])

    intervention_efficacy <- data.frame(
      log_coef = stats::coef(intervention_cox),
      log_coef_se = summary(intervention_cox)$coefficients[, "se(coef)"]
    )

    intervention_efficacy_wtd <- data.frame(
      log_coef = stats::coef(intervention_cox_wtd),
      log_coef_se = summary(intervention_cox_wtd)$coefficients[, "se(coef)"]
    )

  } else if (dtype == "OR") {

    # comparator study
    comparator_efficacy <- agds_eff %>%
      dplyr::filter(rlang::eval_tidy(comparator, .) == comparator.study &
                      rlang::eval_tidy(agds.param.var, .) == agds.param) %>%
    dplyr::transmute(
      log_coef = log(!!agds.estimate),
      log_coef_se = (log(!!agds.ci.upper)
                     - log(!!agds.ci.lower))
      / (2 * stats::qnorm(1 - (1 - CIw) / 2))
    )

    # intervention study to obtain the OR (Odds Ratio)
    formula <- stats::as.formula(
      paste(rlang::as_string(rlang::get_expr(response)), " ~ ",
            rlang::as_string(rlang::get_expr(intervention.arm)), " + ",
            if (!is.null(stralist)) paste(unlist(strsplit(stralist, ", ")),
                                          collapse = " + ") else ""))

    intervention_efficacy <-
      stats::glm(
        formula,
        data = ipds_wts,
        subset = if (!is.null(ipds.param) && !is.null(ipds.param.var)) {
          rlang::eval_tidy(rlang::enquo(ipds.param.var),
                           ipds_wts) == ipds.param
        } else {
          TRUE
        },
        family = stats::binomial(link = "logit")) %>%
      broom::tidy() %>%
      dplyr::filter(
        stringr::str_detect(
          term,
          paste0("^",
                 rlang::as_string(rlang::get_expr(intervention.arm))))) %>%
      dplyr::transmute(log_coef = estimate,
                       log_coef_se = std.error)

    intervention_efficacy_wtd <- suppressWarnings(
      stats::glm(
        formula,
        data = ipds_wts,
        subset = if (!is.null(ipds.param) && !is.null(ipds.param.var)) {
          rlang::eval_tidy(rlang::enquo(ipds.param.var),
                           ipds_wts) == ipds.param
        } else {
          TRUE
        },
        family = stats::binomial(link = "logit"),
        weights = ipds_wts[[rlang::as_name(wt.col)]]) %>%
        broom::tidy() %>%
        dplyr::filter(
          stringr::str_detect(
            term,
            paste0("^",
                   rlang::as_string(rlang::get_expr(intervention.arm))))) %>%
        dplyr::transmute(log_coef = estimate,
                         log_coef_se = std.error))

  } else {
    message("The value of dtype should be HR or OR")
  }

  itc_unadjusted <- estimate_itc(
    intervention = intervention_efficacy,
    comparator = comparator_efficacy,
    outcome = "unadjusted"
  )

  itc_weighted <- estimate_itc(
    intervention = intervention_efficacy_wtd,
    comparator = comparator_efficacy,
    outcome = "weighted"
  )

  itc_result <- dplyr::bind_rows(itc_unadjusted, itc_weighted)

  return(itc_result)
}
