#' @title Summarize the Distribution of Weight Values
#'
#' @param ipds_wts A data frame containing individual patient data from the
#' intervention study, with a column containing the estimated weights (derived
#' using \code{\link{estimate_weights}}).
#' @param intervention.arm The name of the grouping column in the data frame
#' specified by *ipds*, e.g., intervention.arm = TRT. The default is TRT.
#' @param wt.col The name of the estimated weights column in the data frame
#' specified by *ipds_wts*. The default is wt.
#' @param rswt.col The name of the estimated rescaled weights column in the
#' data frame specified by *ipds_wts*. The default is wt_rs.
#' @param digits Specify the number of decimal places for the output results.
#'
#' @importFrom rlang enquo as_name
#' @importFrom assertthat assert_that
#' @importFrom dplyr group_by summarise mutate bind_rows arrange
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom stats median sd quantile
#'
#' @return A data frame containing a summary table of weights and rescaled
#' weights.
#' @export
#'
#' @examples
#' \donttest{
#' summarize_weights(ipds_wts = pts, intervention.arm = TRT)
#' }
#'
#' @name summarize_weights
NULL

# 声明全局变量
utils::globalVariables(c("TRT", "wt", "wt_rs", "stat", "value"))

summarize_weights <- function(ipds_wts, intervention.arm = TRT,
                              wt.col = wt, rswt.col = wt_rs, digits = 2) {

  # 使用enquo捕获参数符号
  intervention.arm <- rlang::enquo(intervention.arm)
  wt.col <- rlang::enquo(wt.col)
  rswt.col <- rlang::enquo(rswt.col)

  # checking--------------------------------------------------------------------
  # 检查ipds_wts是否为data.frame类型
  assertthat::assert_that(
    is.data.frame(ipds_wts),
    msg = "'ipds_wts' is expected to be a data frame with weights")

  # 检查intervention.arm/wt/rswt是否在ipds_wts中
  assertthat::assert_that(
    rlang::as_name(intervention.arm) %in% names(ipds_wts),
    msg = paste0(rlang::as_name(intervention.arm),
                 " can not be found in ipds_wts"))
  assertthat::assert_that(
    rlang::as_name(wt.col) %in% names(ipds_wts),
    msg = paste0(rlang::as_name(wt.col),
                 " can not be found in ipds_wts"))
  assertthat::assert_that(
    rlang::as_name(rswt.col) %in% names(ipds_wts),
    msg = paste0(rlang::as_name(rswt.col),
                 " can not be found in ipds_wts"))

  format_num <- function(x) sprintf(paste0("%.", digits, "f"), x)

  # 计算weights的统计量
  wt_stats <- ipds_wts %>%
    dplyr::group_by(!!intervention.arm) %>%
    dplyr::summarise(
      mean = mean(!!wt.col, na.rm = TRUE) %>% format_num(),
      sd = stats::sd(!!wt.col, na.rm = TRUE) %>% format_num(),
      median = stats::median(!!wt.col, na.rm = TRUE) %>% format_num(),
      Q1 = stats::quantile(!!wt.col, 0.25, na.rm = TRUE) %>% format_num(),
      Q3 = stats::quantile(!!wt.col, 0.75, na.rm = TRUE) %>% format_num(),
      min = min(!!wt.col, na.rm = TRUE) %>% format_num(),
      max = max(!!wt.col, na.rm = TRUE) %>% format_num()
    ) %>%
    tidyr::pivot_longer(cols = -!!intervention.arm,
                        names_to = "stat",
                        values_to = "value") %>%
    dplyr::mutate(var = "Weights")

  # 计算rescaled weights的统计量
  rswt_stats <- ipds_wts %>%
    dplyr::group_by(!!intervention.arm) %>%
    dplyr::summarise(
      mean = mean(!!rswt.col, na.rm = TRUE) %>% format_num(),
      sd = stats::sd(!!rswt.col, na.rm = TRUE) %>% format_num(),
      median = stats::median(!!rswt.col, na.rm = TRUE) %>% format_num(),
      Q1 = stats::quantile(!!rswt.col, 0.25, na.rm = TRUE) %>% format_num(),
      Q3 = stats::quantile(!!rswt.col, 0.75, na.rm = TRUE) %>% format_num(),
      min = min(!!rswt.col, na.rm = TRUE) %>% format_num(),
      max = max(!!rswt.col, na.rm = TRUE) %>% format_num()
    ) %>%
    tidyr::pivot_longer(cols = -!!intervention.arm,
                        names_to = "stat",
                        values_to = "value") %>%
    dplyr::mutate(var = "Rescaled weights")

  stats <- dplyr::bind_rows(wt_stats, rswt_stats) %>%
    tidyr::pivot_wider(names_from = stat, values_from = value) %>%
    dplyr::arrange(!!intervention.arm)

  return(stats)
}
