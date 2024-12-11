#' @title Check Whether the Variables are Balanced After Weighting
#'
#' @param ipds_wts A data frame containing individual patient data from the
#' intervention study, with a column containing the estimated weights (derived
#' using \code{\link{estimate_weights}}).
#' @param agds A data frame containing aggregate summary data from the
#' comparator study.
#' @param summary.list A character list with two elements giving the names of
#' variables for summarizing: the first is a vector of binary variables, and
#' the second is a vector of continuous variables. The variable names must
#' match the column names in *ipds* and do not need to be the same as those
#' in *matching.list*. Use c() if a type is absent.
#' @param matching.list A character list with two elements giving the names of
#' variables for matching: the first is a vector of binary variables, and the
#' second is a vector of continuous variables. The variable names must match
#' the column names in *ipds* and *agds*. Use c() if a type is absent.
#' @param intervention.arm The name of the grouping column in the data frame
#' specified by *ipds*, e.g., `intervention.arm = TRT`. The default is `TRT`.
#' @param comparator The name of the study column in the data frame specified
#' by *agds*, e.g., `comparator = STUDY`. The default is `STUDY`.
#' @param comparator.study A character specifying the comparator study, which
#' must be quoted and cannot be empty (e.g., `comparator.study = "Study XX-1"`).
#' This is the value of the study column in *agds* set by the *comparator*
#' parameter.
#' @param comparator.arm The name of the grouping column in the data frame
#' specified by *agds*, e.g., `comparator.arm = TRT`. The default is `TRT`.
#' @param comparator.n The name of the subjects number column in the data frame
#' specified by *agds*, e.g., `comparator.n = N`. The default is `N`.
#' @param wt.col The name of the estimated weights column in the data frame
#' specified by *ipds_wts*. The default is `wt`.
#'
#' @importFrom rlang enquo as_name
#' @importFrom assertthat assert_that
#' @importFrom dplyr union filter mutate select summarise rename group_by
#' @importFrom dplyr distinct bind_rows left_join transmute matches across
#' @importFrom magrittr %>%
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect starts_with all_of
#' @importFrom data.table :=
#'
#' @return A data frame containing all specified variables summarised before
#' and after weighting.
#' @export
#'
#' @examples
#' \donttest{
#' cov <- list(
#'   binary = c("ECOG", "SMK", "METBRAIN"),
#'   continuous = c("BMI", "DIAG")
#' )
#'
#' cov_all <- list(
#'   binary = c("SEX", "ECOG", "SMK", "METBRAIN", "METLIVER"),
#'   continuous = c("BMI", "DIAG", "WEIGHT", "HEIGHT")
#' )
#'
#' baseline <- check_matching(
#'   ipds_wts = pts, agds = AgD_bl,
#'   summary.list = cov_all, matching.list = cov,
#'   intervention.arm = TRT,
#'   comparator = STUDY, comparator.study = "Study XX-1",
#'   comparator.n = N, comparator.arm = TRT)
#'
#' baseline
#' }
#'
#' @name check_matching
NULL

# 声明全局变量
utils::globalVariables(c("TRT", "STUDY", "N", "wt", "intervention_n",
                         "intervention_ess", "src", "total", "p", "n",
                         "type", "variable", "statistic", "value",
                         ".data"))

check_matching <- function(ipds_wts, agds,
                           summary.list, matching.list = summary.list,
                           intervention.arm = TRT,
                           comparator = STUDY, comparator.study,
                           comparator.arm = TRT, comparator.n = N,
                           wt.col = wt) {

  # 使用enquo捕获参数符号
  intervention.arm <- rlang::enquo(intervention.arm)
  comparator <- rlang::enquo(comparator)
  comparator.arm <- rlang::enquo(comparator.arm)
  comparator.n <- rlang::enquo(comparator.n)
  wt.col <- rlang::enquo(wt.col)

  # checking--------------------------------------------------------------------
  # 检查ipds和agds是否为data.frame类型
  assertthat::assert_that(
    is.data.frame(ipds_wts),
    msg = "'ipds_wts' is expected to be a data frame with wt")
  assertthat::assert_that(
    is.data.frame(agds),
    msg = "'agds' is expected to be a data frame")

  # 检查intervention.arm和comparator.arm是否在ipds和agds中
  assertthat::assert_that(
    rlang::as_name(intervention.arm) %in% names(ipds_wts),
    msg = paste0(rlang::as_name(intervention.arm),
                 " can not be found in ipds_wts"))
  assertthat::assert_that(
    rlang::as_name(comparator.arm) %in% names(agds),
    msg = paste0(rlang::as_name(comparator.arm),
                 " can not be found in agds"))

  # 检查matching.list的类型及其变量是否在数据框中
  summary_cov <- dplyr::union(unlist(summary.list[[1]]),
                              unlist(summary.list[[2]]))

  assertthat::assert_that(
    length(summary_cov) != 0,
    msg = "summary.list is NULL")

  assertthat::assert_that(
    is.character(summary_cov),
    msg = "matching variables are expected to be a character vector")
  missing_vars <- summary_cov[!(summary_cov %in% names(ipds_wts))]
  if (length(missing_vars) > 0) {
    stop(paste(missing_vars, collapse = ", "), " cannot be found in ipds_wts")
  }

  # 查找连续变量的mean是否在 agds 中
  # if (length(summary.list[[2]]) > 0) {
  #   summary_cov_cont_mean <- paste0(unlist(summary.list[[2]]), ".mean")
  #   missing_vars_cont_mean <-
  #     summary_cov_cont_mean[!summary_cov_cont_mean %in% names(agds)]
  #   if (length(missing_vars_cont_mean) > 0) {
  #     stop(paste(missing_vars_cont_mean, collapse = "; "),
  #          " cannot be found in agds")
  #   }
  # }

  # 保存intervention.arm类型和顺序 ---------------------------------------------
  intervention_levels <-
    levels(ipds_wts[[rlang::as_name(intervention.arm)]])
  intervention_arm_as_fac <-
    is.factor(ipds_wts[[rlang::as_name(intervention.arm)]])

  match_cov <- dplyr::union(unlist(matching.list[[1]]),
                            unlist(matching.list[[2]]))

  summary_cov_binary <- unlist(summary.list[[1]])
  summary_cov_binary_comp <-
    summary_cov_binary[summary_cov_binary %in% colnames(agds)]

  summary_cov_continuous <- unlist(summary.list[[2]])

  ess <- estimate_ess(
    ipds_wts = ipds_wts, agds = agds,
    intervention.arm = !!intervention.arm,
    comparator = !!comparator,
    comparator.study = comparator.study,
    comparator.arm = !!comparator.arm,
    comparator.n = !!comparator.n) %>%
    dplyr::filter(!!intervention.arm != "total") %>%
    dplyr::select(!!intervention.arm,
                  comparator,
                  intervention_pre = intervention_n,
                  intervention_post = intervention_ess) %>%
    tidyr::pivot_longer(
      cols = tidyselect::starts_with("comparator") |
        tidyselect::starts_with("intervention"),
      names_to = c("src"),
      values_to = "total") %>%
    dplyr::select(!!intervention.arm, src, total)


  # 二分类终点
  cal_mean_binary <- function(df, vars, weight = NULL) {
    if (!is.null(weight)) {
      df %>% dplyr::summarise(
        dplyr::across(
          tidyselect::all_of(vars),
          ~ weighted.mean(., get(weight), na.rm = TRUE)))
    } else {
      df %>% dplyr::summarise(
        dplyr::across(
          tidyselect::all_of(vars),
          ~ mean(., na.rm = TRUE)))
    }
  }

  comp_binary <- agds %>%
    dplyr::rename(!!rlang::as_name(intervention.arm)
                  := !!rlang::as_name(comparator.arm)) %>%
    dplyr::filter(!!comparator == comparator.study) %>%
    dplyr::group_by(!!intervention.arm) %>%
    cal_mean_binary(summary_cov_binary_comp) %>%
    dplyr::distinct() %>%
    tidyr::pivot_longer(
      cols = tidyselect::all_of(summary_cov_binary_comp),
      names_to = "variable",
      values_to = "p") %>%
    dplyr::mutate(src = "comparator")

  # 计算干预组未加权的数据均值
  inter_pre_binary <- ipds_wts %>%
    dplyr::group_by(!!intervention.arm) %>%
    cal_mean_binary(summary_cov_binary) %>%
    dplyr::distinct() %>%
    tidyr::pivot_longer(
      cols = tidyselect::all_of(summary_cov_binary),
      names_to = "variable", values_to = "p") %>%
    dplyr::mutate(src = "intervention_pre")

  # 计算干预组加权的数据均值
  inter_post_binary <- ipds_wts %>%
    dplyr::group_by(!!intervention.arm) %>%
    cal_mean_binary(summary_cov_binary, weight = rlang::as_name(wt.col)) %>%
    dplyr::distinct() %>%
    tidyr::pivot_longer(
      cols = tidyselect::all_of(summary_cov_binary),
      names_to = "variable",
      values_to = "p") %>%
    dplyr::mutate(src = "intervention_post")

  combine_binary <- dplyr::bind_rows(
    inter_pre_binary,
    inter_post_binary,
    comp_binary) %>%
    dplyr::left_join(ess, by = c(rlang::as_name(intervention.arm), "src")) %>%
    dplyr::mutate(type = "binary",
                  n = round(p * total, digits = 0),
                  p_oth = 1 - p,
                  n_oth = total - n) %>%
    tidyr::pivot_longer(
      cols = c("n", "p", "n_oth", "p_oth"),
      names_to = "statistic",
      values_to = "value") %>%
    dplyr::transmute(type, !!intervention.arm, variable,
                     adj = ifelse(variable %in% match_cov, "Y", "N"),
                     src, total, statistic, value)


  # 连续型终点
  cal_mean_continuous <- function(df, vars, weight = NULL) {
    if (!is.null(weight)) {
      df %>% dplyr::summarise(dplyr::across(tidyselect::all_of(vars), list(
        mean = ~ weighted.mean(.x, get(weight), na.rm = TRUE),
        sd = ~ sqrt(weighted.mean((.x - weighted.mean(.x, get(weight),
                                                      na.rm = TRUE))^2,
                                  get(weight),
                                  na.rm = TRUE))
      )))
    } else {
      df %>% dplyr::summarise(dplyr::across(tidyselect::all_of(vars), list(
        mean = ~ mean(.x, na.rm = TRUE),
        sd = ~ sd(.x, na.rm = TRUE))))
    }
  }

  comp_continuous <- agds %>%
    dplyr::rename(!!rlang::as_name(intervention.arm)
                  := !!rlang::as_name(comparator.arm)) %>%
    dplyr::filter(!!comparator == comparator.study) %>%
    dplyr::select(!!intervention.arm, dplyr::matches("\\.(mean|sd)$")) %>%
    tidyr::pivot_longer(
      cols = -!!intervention.arm,
      names_to = c("variable", "statistic"),
      names_sep = "\\.",
      values_to = "value"
    ) %>%
    dplyr::mutate(src = "comparator")

  # 计算干预组未加权的数据均值
  inter_pre_continuous <- ipds_wts %>%
    dplyr::group_by(!!intervention.arm) %>%
    cal_mean_continuous(summary_cov_continuous) %>%
    dplyr::distinct() %>%
    tidyr::pivot_longer(
      cols = -!!intervention.arm,
      names_to = c("variable", "statistic"),
      names_pattern = "(.*)_(.*)",
      values_to = "value"
    ) %>%
    dplyr::mutate(src = "intervention_pre")

  # 计算干预组加权的数据均值
  inter_post_continuous <- ipds_wts %>%
    dplyr::group_by(!!intervention.arm) %>%
    cal_mean_continuous(summary_cov_continuous,
                        weight = rlang::as_name(wt.col)) %>%
    dplyr::distinct() %>%
    tidyr::pivot_longer(
      cols = -!!intervention.arm,
      names_to = c("variable", "statistic"),
      names_pattern = "(.*)_(.*)",
      values_to = "value"
    ) %>%
    dplyr::mutate(src = "intervention_post")

  combine_continuous <- dplyr::bind_rows(
    inter_pre_continuous,
    inter_post_continuous,
    comp_continuous) %>%
    dplyr::left_join(ess, by = c(rlang::as_name(intervention.arm), "src")) %>%
    dplyr::transmute(type = "continuous", !!intervention.arm, variable,
                     adj = ifelse(variable %in% match_cov, "Y", "N"),
                     src, total, statistic, value)

  combine <- dplyr::bind_rows(combine_binary, combine_continuous)

  # 恢复intervention.arm的原始类型和顺序
  if (intervention_arm_as_fac) {
    combine <- combine %>%
      mutate(!!rlang::as_name(intervention.arm)
             := factor(.data[[rlang::as_name(intervention.arm)]],
                       levels = intervention_levels)) %>%
      arrange(!!intervention.arm)
  }

  return(combine)
}
