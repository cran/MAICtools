#' @title Conduct non-Anchored Matching-Adjusted Indirect Comparison (MAIC).
#'
#' @param unds_wts A combined data frame containing individual efficacy data
#' from the intervention study and pseudo efficacy data from the comparator
#' study.
#' @param unds.arm The name of the grouping column in the combined data frame
#' specified by *unds_wts*, e.g., comparator.arm = TRT. The default is TRT.
#' @param comparator.study A character specifying or presenting the comparator
#' study, e.g., comparator.study = "Study XX-1".
#' @param unds.param.var The name of the column that specifies only a subset
#' of the rows of the data to be used.
#' @param unds.param A character specifying the subset of the rows to be used.
#' This is the value of the column set by the *unds.param.var*.
#' @param time The name of the survival or follow up time column in the
#' combined data frame.
#' @param status The status indicator, normally 0 = event, 1 = censored.
#' Can be reseted using the *event* parameter.
#' @param event A numeric value that represents the survival status, 0 = event,
#' 1 = censored.
#' @param response The name of the response status column in the *unds_wts*.
#' @param dtype Two options are available: "HR" and "OR". The default is "HR".
#' @param wt.col The name of the estimated weights column in the data frame
#' specified by *unds_wts*. The default is wt.
#' @param CIw The numeric value specifying the width of the confidence interval,
#' with a default of 0.95.
#' @param digits Specify the number of decimal places for the output results.
#'
#' @importFrom rlang enquo as_name eval_tidy get_expr as_string
#' @importFrom assertthat assert_that
#' @importFrom dplyr transmute bind_rows filter
#' @importFrom stats glm binomial as.formula confint
#' @importFrom survival Surv coxph
#' @importFrom broom tidy
#' @importFrom stringr str_detect
#' @importFrom tibble as_tibble
#' @importFrom data.table :=
#'
#' @return A data frame containing the non-anchored matching-adjusted indirect
#' comparison results.
#' @export
#'
#' @examples
#' \donttest{
#' results3 <- unanchored_maic(
#'   unds_wts = unpts, unds.arm = ARM,
#'   comparator.study = "Study XX-1",
#'   unds.param.var = PARAMCD, unds.param = "OS",
#'   time = AVAL, status = CNSR, event = 0,
#'   dtype = "HR")
#'
#' results3
#'
#' results4 <- unanchored_maic(
#'   unds_wts = unpts, unds.arm = ARM,
#'   #' unds.param = "ORR",
#'   #' comparator.study = "Study XX-1",
#'   response = CNSR,
#'   dtype = "OR")
#'
#' results4
#' }
#'
#' @name unanchored_maic
NULL

# 声明全局变量
utils::globalVariables(c("ARM", "wt", "exp(coef)", "term", "estimate",
                         "conf.low", "conf.high"))

unanchored_maic <- function(unds_wts, unds.arm = ARM,
                            comparator.study = NULL,
                            unds.param.var = NULL, unds.param = NULL,
                            time = NULL, status = NULL, event = 0,
                            response = NULL,
                            dtype = "HR", wt.col = wt, CIw = 0.95,
                            digits = 2) {

  # 使用 enquo 捕获参数符号
  unds.arm <- rlang::enquo(unds.arm)
  unds.param.var <- rlang::enquo(unds.param.var)
  time <- rlang::enquo(time)
  status <- rlang::enquo(status)
  response <- rlang::enquo(response)
  wt.col <- rlang::enquo(wt.col)

  # checking--------------------------------------------------------------------
  # 检查unds_wts是否为data.frame类型
  assertthat::assert_that(
    is.data.frame(unds_wts),
    msg = "'unds_wts' is expected to be a data frame")

  # 检查变量unds.arm、wt.col是否在数据集unds_wts中
  assertthat::assert_that(
    rlang::as_name(unds.arm) %in% names(unds_wts),
                          msg = paste0(rlang::as_name(unds.arm),
                                       " can not be found in unds_wts"))
  assertthat::assert_that(
    rlang::as_name(wt.col) %in% names(unds_wts),
                          msg = paste0(rlang::as_name(wt.col),
                                       " can not be found in unds_wts"))

  # 当dtype = "HR"时，检查变量time、status是否在数据集unds_wts中
  if (dtype == "HR") {
    assertthat::assert_that(
      rlang::as_name(time) %in% names(unds_wts),
                            msg = paste0(rlang::as_name(time),
                                         " can not be found unds_wts"))
    assertthat::assert_that(
      rlang::as_name(status) %in% names(unds_wts),
                            msg = paste0(rlang::as_name(status),
                                         " can not be found in unds_wts"))
  }

  # 当dtype = "OR"时，检查变量response是否在数据集unds_wts中
  if (dtype == "OR") {
    assertthat::assert_that(
      rlang::as_name(response) %in% names(unds_wts),
                            msg = paste0(rlang::as_name(response),
                                         " can not be found in unds_wts"))
  }


  # ----------------------------------------------------------------------------
  lower_ci <- paste0("lower_", CIw, "ci")
  upper_ci <- paste0("upper_", CIw, "ci")

  if (dtype == "HR") {
    # obtain the HR (Hazard Ratio)
    formula <- stats::as.formula(
      paste0("survival::Surv(", rlang::as_string(rlang::get_expr(time)), ", ",
             rlang::as_string(rlang::get_expr(status)), " == ", event, ") ~ ",
             rlang::as_string(rlang::get_expr(unds.arm))))

    unadjusted_cox <- survival::coxph(
      formula, data = unds_wts,
      subset = if (!is.null(unds.param) && !is.null(unds.param.var)) {
        rlang::eval_tidy(rlang::enquo(unds.param.var), unds_wts) == unds.param
      } else {
        TRUE
      },
      method = "efron")

    weighted_cox <- survival::coxph(
      formula, data = unds_wts,
      subset = if (!is.null(unds.param) && !is.null(unds.param.var)) {
        rlang::eval_tidy(rlang::enquo(unds.param.var), unds_wts) == unds.param
      } else {
        TRUE
      },
      method = "efron",
      weights = unds_wts[[rlang::as_name(wt.col)]])

    unadjusted_efficacy <- summary(unadjusted_cox)$conf.int %>%
      as.data.frame() %>%
      dplyr::transmute(
        comparator = comparator.study,
        param = unds.param,
        type = "un-anchored",
        outcome = "unadjusted",
        effects = sprintf(paste0("%.", digits, "f"), `exp(coef)`),
        !!lower_ci := sprintf(paste0("%.", digits, "f"),
                              exp(stats::confint(unadjusted_cox,
                                                 level = CIw))[1]),
        !!upper_ci := sprintf(paste0("%.", digits, "f"),
                              exp(stats::confint(unadjusted_cox,
                                                 level = CIw))[2])) %>%
      tibble::as_tibble(rownames = NULL)

    weighted_efficacy <- summary(weighted_cox)$conf.int %>%
      as.data.frame() %>%
      dplyr::transmute(
        comparator = comparator.study,
        param = unds.param,
        type = "un-anchored",
        outcome = "weighted",
        effects = sprintf(paste0("%.", digits, "f"), `exp(coef)`),
        !!lower_ci := sprintf(paste0("%.", digits, "f"),
                              exp(stats::confint(weighted_cox,
                                                 level = CIw))[1]),
        !!upper_ci := sprintf(paste0("%.", digits, "f"),
                              exp(stats::confint(weighted_cox,
                                                 level = CIw))[2])) %>%
      tibble::as_tibble(rownames = NULL)

    itc_result <- dplyr::bind_rows(unadjusted_efficacy, weighted_efficacy)

  } else if (dtype == "OR") {
    # obtain the OR (Odds Ratio)
    formula <- stats::as.formula(
      paste0(rlang::as_string(rlang::get_expr(response)), " ~ ",
             rlang::as_string(rlang::get_expr(unds.arm))))

    unadjusted_glm <-
      stats::glm(
        formula,
        data = unds_wts,
        subset = if (!is.null(unds.param) && !is.null(unds.param.var)) {
          rlang::eval_tidy(rlang::enquo(unds.param.var),
                           unds_wts) == unds.param
        } else {
          TRUE
        },
        family = stats::binomial(link = "logit"))

    weighted_glm <- suppressWarnings(
      stats::glm(
        formula,
        data = unds_wts,
        subset = if (!is.null(unds.param) && !is.null(unds.param.var)) {
          rlang::eval_tidy(rlang::enquo(unds.param.var),
                           unds_wts) == unds.param
        } else {
          TRUE
        },
        family = stats::binomial(link = "logit"),
        weights = unds_wts[[rlang::as_name(wt.col)]]))

    unadjusted_efficacy <- suppressWarnings(
      broom::tidy(unadjusted_glm, conf.int = TRUE, conf.level = CIw) %>%
        dplyr::filter(stringr::str_detect(
          term, paste0("^",
                       rlang::as_string(rlang::get_expr(unds.arm))))) %>%
        dplyr::transmute(
          comparator = comparator.study,
          param = unds.param,
          type = "un-anchored",
          outcome = "unadjusted",
          effects = sprintf(paste0("%.", digits, "f"), exp(estimate)),
          !!lower_ci := sprintf(paste0("%.", digits, "f"), exp(conf.low)),
          !!upper_ci := sprintf(paste0("%.", digits, "f"), exp(conf.high))))

    weighted_efficacy <- suppressWarnings(
      broom::tidy(weighted_glm, conf.int = TRUE, conf.level = CIw) %>%
        dplyr::filter(stringr::str_detect(
          term, paste0("^", rlang::as_string(rlang::get_expr(unds.arm))))) %>%
        dplyr::transmute(
          comparator = comparator.study,
          param = unds.param,
          type = "un-anchored",
          outcome = "weighted",
          effects = sprintf(paste0("%.", digits, "f"), exp(estimate)),
          !!lower_ci := sprintf(paste0("%.", digits, "f"), exp(conf.low)),
          !!upper_ci := sprintf(paste0("%.", digits, "f"), exp(conf.high))))

    itc_result <- dplyr::bind_rows(unadjusted_efficacy, weighted_efficacy)
  } else {
    message("The value of dtype should be HR or OR")
  }

  return(itc_result)
}

