#' @title Estimate Effective Sample Size (ESS)
#'
#' @param ipds_wts A data frame containing individual patient data from the
#' intervention study, with a column containing the estimated weights (derived
#' using \code{\link{estimate_weights}}).
#' @param agds A data frame containing aggregate summary data from the
#' comparator study.
#' @param intervention.arm The name of the grouping column in the data frame
#' specified by *ipds*, e.g., intervention.arm = TRT. The default is TRT.
#' @param comparator The name of the study column in the data frame specified
#' by *agds*, e.g., comparator = STUDY. The default is STUDY.
#' @param comparator.study A character specifying the comparator study, which
#' must be quoted and cannot be empty (e.g., comparator.study = "Study XX-1").
#' This is the value of the study column in *agds* set by the *comparator*
#' parameter.
#' @param comparator.arm The name of the grouping column in the data frame
#' specified by *agds*, e.g., comparator.arm = TRT. The default is TRT.
#' @param comparator.n A The name of the subjects number column in the data
#' frame specified by *agds*, e.g., comparator.n = N. The default is N.
#' @param wt.col The name of the estimated weights column in the data frame
#' specified by *ipds_wts*. The default is wt.
#' @param digits Specify the number of decimal places for the output results.
#'
#' @importFrom rlang enquo as_name
#' @importFrom assertthat assert_that
#' @importFrom dplyr rename filter bind_rows pull
#' @importFrom magrittr %>%
#' @importFrom purrr map
#' @importFrom tibble tibble
#' @importFrom data.table :=
#'
#' @return A data frame containing effective sample size (ESS) after weighting.
#' @export
#'
#' @examples
#' \donttest{
#' ess <- estimate_ess(
#'   ipds_wts = pts, agds = AgD_bl,
#'   intervention.arm = TRT,
#'   comparator = STUDY, comparator.study = "Study XX-1", comparator.arm = TRT,
#'   comparator.n = N)
#' ess
#' }
#'
#' @name estimate_ess

# 声明全局变量
utils::globalVariables(c("TRT", "STUDY", "N", "wt", ".data"))

estimate_ess <- function(ipds_wts, agds,
                         intervention.arm = TRT,
                         comparator = STUDY,
                         comparator.study = NULL,
                         comparator.arm = TRT,
                         comparator.n = N,
                         wt.col = wt,
                         digits = 1) {

  # 将参数捕获为quosure
  intervention.arm <- rlang::enquo(intervention.arm)
  comparator <- rlang::enquo(comparator)
  comparator.arm <- rlang::enquo(comparator.arm)
  comparator.n <- rlang::enquo(comparator.n)
  wt.col <- rlang::enquo(wt.col)

  # 检查输入的数据集类型
  assertthat::assert_that(
    is.data.frame(ipds_wts),
    msg = "'ipds_wts' is expected to be a data frame with weights")
  assertthat::assert_that(
    is.data.frame(agds),
    msg = "'agds' is expected to be a data frame")

  # 检查变量是否存在于对应的数据集中
  assertthat::assert_that(
    rlang::as_name(intervention.arm) %in% names(ipds_wts),
    msg = paste0(rlang::as_name(intervention.arm),
                 " can not be found in ipds_wts"))
  assertthat::assert_that(
    rlang::as_name(wt.col) %in% names(ipds_wts),
    msg = paste0(rlang::as_name(wt.col),
                 " can not be found in ipds_wts"))
  assertthat::assert_that(
    rlang::as_name(comparator.arm) %in% names(agds),
    msg = paste0(rlang::as_name(comparator.arm),
                 " can not be found in agds"))
  assertthat::assert_that(
    rlang::as_name(comparator.n) %in% names(agds),
    msg = paste0(rlang::as_name(comparator.n),
                 " can not be found in agds"))

  # 保存intervention.arm的原始类型
  intervention_levels <-
    c(levels(ipds_wts[[rlang::as_name(intervention.arm)]]), "total")
  intervention_arm_as_fac <-
    is.factor(ipds_wts[[rlang::as_name(intervention.arm)]])

  # 重命名agds的变量
  agds <- dplyr::rename(agds,
                        !!rlang::as_name(intervention.arm)
                        := !!rlang::as_name(comparator.arm))

  # 获取唯一的组名
  groups <- unique(c(ipds_wts[[rlang::as_name(intervention.arm)]],
                     agds[[rlang::as_name(intervention.arm)]]))

  results_grp <- purrr::map(groups, function(group) {
    intervention_grp <- dplyr::filter(ipds_wts, !!intervention.arm == group)
    comparator_n <- agds %>% dplyr::filter(
      !!comparator == comparator.study & !!intervention.arm == group) %>%
      dplyr::pull(!!comparator.n)

    if (is.null(comparator_n)
        || (length(comparator_n) == 1 && is.na(comparator_n))) {
      return(NULL)
    }

    intervention_n <- nrow(intervention_grp)
    intervention_ess <- round(
      sum(intervention_grp[[rlang::as_name(wt.col)]])^2
      / sum(intervention_grp[[rlang::as_name(wt.col)]]^2),
      digits = 0)
    ratio_ess <- sprintf(paste0("%.", digits, "f"),
                         intervention_ess / intervention_n * 100)

    tibble::tibble(
      !!rlang::as_name(intervention.arm) := group,
      comparator = comparator_n,
      intervention_n = intervention_n,
      intervention_ess = intervention_ess,
      ratio_ess = ratio_ess
    )
  }) %>% dplyr::bind_rows()

  # 不分组的结果
  intervention_ugrp <- ipds_wts
  comparator_n_ugrp <- agds %>%
    dplyr::filter(!!comparator == comparator.study) %>%
    dplyr::pull(!!comparator.n) %>% sum(na.rm = TRUE)
  intervention_n_ugrp <- nrow(intervention_ugrp)
  intervention_ess_ugrp <-
    round(sum(intervention_ugrp[[rlang::as_name(wt.col)]])^2
          / sum(intervention_ugrp[[rlang::as_name(wt.col)]]^2), digits = 0)
  ratio_ess_ugrp <-
    sprintf(paste0("%.", digits, "f"),
            intervention_ess_ugrp / intervention_n_ugrp * 100)

  results_ugrp <- tibble::tibble(
    !!rlang::as_name(intervention.arm) := "total",
    comparator = comparator_n_ugrp,
    intervention_n = intervention_n_ugrp,
    intervention_ess = intervention_ess_ugrp,
    ratio_ess = ratio_ess_ugrp
  )

  results <- dplyr::bind_rows(results_grp, results_ugrp)

  # 恢复intervention.arm的原始类型和顺序
  if (intervention_arm_as_fac) {
    results <- results %>%
      mutate(!!rlang::as_name(intervention.arm)
             := factor(.data[[rlang::as_name(intervention.arm)]],
                       levels = intervention_levels)) %>%
      arrange(!!intervention.arm)
  }

  return(results)
}
