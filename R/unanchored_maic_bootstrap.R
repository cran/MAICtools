#' @title Conduct non-Anchored Matching-Adjusted Indirect Comparison (MAIC) and
#' Calculate Confidence Intervals (CIs) Using Bootstrap.
#'
#' @description
#'     Two different methods for estimating a 95% confidence interval (CI) from
#'     the bootstrap samples were explored:
#'     * Percentile CIs
#'     * Bias-corrected and accelerated (BCa) CIs
#'
#' @param ipds A data frame containing individual patient data from the
#' intervention study, with baseline characteristic variables for matching.
#' @param psds A data frame containing pseudo data from the comparator study.
#' @param agds A data frame containing aggregate summary data from the
#' comparator study.
#' @param matching.list A character list with two elements giving the names
#' of variables for matching: the first is a vector of binary variables, and
#' the second is a vector of continuous variables. The variable names must
#' match the column names in *ipds* and *agds*. Use c() if a type is absent.
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
#' @param ipds.param.var The name of the column that specifies only a subset
#' of the *ipds* to be used.
#' @param ipds.param A character specifying the subset of the rows to be used.
#' This is the value of the column set by the *ipds.param.var*.
#' @param psds.param.var The name of the column that specifies only a specifyed
#' result of the *psds* to be used.
#' @param psds.param A character specifying the subset of the rows to be used.
#' This is the value of the column set by the *psds.param.var*.
#' @param time The name of the survival or follow up time column.
#' @param status The status indicator, normally 0 = event, 1 = censored. Can be
#' reseted using the *event* parameter.
#' @param event A numeric value that represents the survival status, 0 = event,
#'  1 = censored.
#' @param response The name of the response status column.
#' @param dtype Two options are available: "HR" and "OR". The default is "HR".
#' @param n.samples The number of bootstrap replicates.
#' @param CIw The numeric value specifying the width of the confidence interval,
#'  with a default of 0.95.
#' @param digits Specify the number of decimal places for the output results.
#' @param ... Refer to \link[boot:boot]{boot} for additional parameters.
#'
#' @importFrom rlang enquo as_name eval_tidy get_expr as_string quo_is_null
#' @importFrom assertthat assert_that
#' @importFrom dplyr union filter
#' @importFrom stats glm binomial coef relevel as.formula quantile
#' @importFrom survival Surv coxph
#' @importFrom boot boot boot.ci
#' @importFrom graphics abline hist
#' @importFrom grDevices recordPlot
#'
#' @return A list containing 2 objects.
#' First, a data frame containing the non-anchored matching-adjusted indirect
#' comparison results.
#' Second, a bootstrapping diagnostics histogram.
#' @export
#'
#' @examples
#' \donttest{
#' cov <- list(
#'   c("ECOG", "SMK", "METBRAIN"),
#'   c("BMI", "DIAG")
#' )
#'
#' results5 <- unanchored_maic_bootstrap(
#'   ipds = IPD,
#'   agds = AgD_bl,
#'   psds = pseudo,
#'   matching.list = cov,
#'   intervention.arm = TRT,
#'   comparator = STUDY,
#'   comparator.study = "Study XX-1",
#'   comparator.arm = TRT,
#'   time = AVAL, status = CNSR, event = 0,
#'   dtype = "HR",
#'   ipds.param.var = PARAMCD, ipds.param = "OS",
#'   psds.param.var = NULL, psds.param = NULL,
#'   n.samples = 1000
#' )
#'
#' results5$results
#' results5$plot
#' }
#'
#' @name unanchored_maic_bootstrap
NULL

# 声明全局变量
utils::globalVariables(c("ARM", "wt"))

unanchored_maic_bootstrap <- function(
    ipds, psds, agds, matching.list, intervention.arm,
    comparator, comparator.study, comparator.arm,
    ipds.param.var = NULL, ipds.param = NULL,
    psds.param.var = NULL, psds.param = NULL,
    time = NULL, status = NULL, event = 0, response = NULL,
    dtype = "HR", n.samples, CIw = 0.95, digits = 2,
    ...) {

  # 捕获列名并转换为符号
  intervention.arm <- rlang::enquo(intervention.arm)
  comparator <- rlang::enquo(comparator)
  comparator.arm <- rlang::enquo(comparator.arm)
  ipds.param <- rlang::enquo(ipds.param)
  ipds.param.var <- rlang::enquo(ipds.param.var)
  psds.param <- rlang::enquo(psds.param)
  psds.param.var <- rlang::enquo(psds.param.var)
  time <- rlang::enquo(time)
  status <- rlang::enquo(status)
  response <- rlang::enquo(response)

  # checking--------------------------------------------------------------------
  # 检查 ipds, psds, agds 是否为 data.frame 类型
  assertthat::assert_that(
    is.data.frame(ipds),
    msg = "'ipds' is expected to be a data frame")
  assertthat::assert_that(
    is.data.frame(psds),
    msg = "'psds' is expected to be a data frame")
  assertthat::assert_that(
    is.data.frame(agds),
    msg = "'agds' is expected to be a data frame")

  # 检查 intervention.arm 和 comparator.arm 是否在 ipds 和 agds 中
  assertthat::assert_that(
    rlang::as_name(intervention.arm) %in% names(ipds),
    msg = paste0(rlang::as_name(intervention.arm),
                 " can not be found in ipds"))
  assertthat::assert_that(
    rlang::as_name(comparator.arm) %in% names(agds),
    msg = paste0(rlang::as_name(comparator.arm),
                 " can not be found in agds"))

  # 检查 matching.list 类型及其变量是否在数据框中
  if (!is.list(matching.list) || length(matching.list) < 2) {
    stop("'matching.list' should be a list with at least two elements")
  }

  if (!is.character(unlist(matching.list[[1]]))
      || !is.character(unlist(matching.list[[2]]))) {
    stop("Elements of 'matching.list' should be character vectors")
  }

  match_cov <- dplyr::union(unlist(matching.list[[1]]),
                            unlist(matching.list[[2]]))

  assertthat::assert_that(
    length(match_cov) != 0,
    msg = "matching.list is NULL")
  assertthat::assert_that(
    is.character(match_cov),
    msg = "matching variables are expected to be a character vector")
  missing_vars <- match_cov[!(match_cov %in% names(ipds))]
  if (length(missing_vars) > 0) {
    stop(paste(missing_vars, collapse = ", "), " cannot be found in ipds")
  }

  # 当 dtype = "HR" 时，检查变量 time, status 是否在数据集 ipds, psds 中
  if (dtype == "HR") {
    assertthat::assert_that(
      rlang::as_name(time) %in% names(ipds),
      msg = paste0(rlang::as_name(time),
                   " can not be found ipds"))
    assertthat::assert_that(
      rlang::as_name(status) %in% names(ipds),
      msg = paste0(rlang::as_name(status),
                   " can not be found in ipds"))

    assertthat::assert_that(
      rlang::as_name(time) %in% names(psds),
      msg = paste0(rlang::as_name(time),
                   " can not be found psds"))
    assertthat::assert_that(
      rlang::as_name(status) %in% names(psds),
      msg = paste0(rlang::as_name(status),
                   " can not be found psds"))
  }

  # 当 dtype = "OR" 时，检查变量 response 是否在数据集 ipds, psds 中
  if (dtype == "OR") {
    assertthat::assert_that(
      rlang::as_name(response) %in% names(ipds),
      msg = paste0(rlang::as_name(response),
                   " can not be found ipds"))

    assertthat::assert_that(
      rlang::as_name(response) %in% names(psds),
      msg = paste0(rlang::as_name(response),
                   " can not be found psds"))
  }

  # ----------------------------------------------------------------------------
  lower_ci <- paste0("lower_", CIw, "ci")
  upper_ci <- paste0("upper_", CIw, "ci")

  if (!rlang::quo_is_null(ipds.param) && !rlang::quo_is_null(ipds.param.var)) {
    ipds <- ipds %>% dplyr::filter(!!ipds.param.var == !!ipds.param)
  }
  if (!rlang::quo_is_null(psds.param) && !rlang::quo_is_null(psds.param.var)) {
    psds <- psds %>% dplyr::filter(!!psds.param.var == !!psds.param)
  }

  # 自举统计函数
  bootstrap_func <- function(data, inds) {
    # 使用抽样索引创建自举样本
    bootstrap_ipds <- data[inds, ]

    # 过滤 active arm 的数据
    ipds_active <- bootstrap_ipds %>% filter(!!intervention.arm == "active")

    # 估算权重
    ipds_wt <- estimate_weights(
      ipds = ipds_active,
      agds = agds,
      matching.list = matching.list,
      intervention.arm = !!intervention.arm,
      comparator = !!comparator,
      comparator.study = comparator.study,
      comparator.arm = !!comparator.arm
    ) %>% mutate(ARM = "Intervention")

    # 给比较组数据权重 1
    ps_wt <- psds %>% mutate(wt = 1, ARM = "Comparator")

    # 合并数据
    combined_data <- bind_rows(ipds_wt, ps_wt) %>%
      mutate(ARM = stats::relevel(as.factor(ARM), ref = "Comparator"),
             wt = pmax(wt, 1e-6))  # 避免权重为 0

    if (dtype == "HR") {
      # 获得 HR (Hazard Ratio)
      formula <- stats::as.formula(
        paste0("survival::Surv(",
               rlang::as_string(rlang::get_expr(time)), ", ",
               rlang::as_string(rlang::get_expr(status)),
               " == ", event, ") ~ ", "ARM"))

      cox_model <- survival::coxph(
        formula, data = combined_data, method = "efron", weights = wt)

      effect <- exp(cox_model$coefficients)

    } else if (dtype == "OR") {
      # 获得 OR (Odds Ratio)
      formula <- stats::as.formula(
        paste0(rlang::as_string(rlang::get_expr(response)), " ~ ", "ARM"))

      weighted_glm <- suppressWarnings(stats::glm(
        formula,
        data = combined_data,
        family = stats::binomial(link = "logit"),
        weights = wt))

      effect <- exp(as.numeric(stats::coef(weighted_glm)[2]))

    } else {
      message("The value of dtype should be HR or OR")
    }

    return(effect)
  }

  # 运行自举
  result_bootstraps <- boot::boot(data = ipds,
                                  statistic = bootstrap_func,
                                  R = n.samples, ...)

  # 保证结果是一个数值向量
  if (!is.null(dim(result_bootstraps$t)) && dim(result_bootstraps$t)[2] > 1) {
    stop("bootstrap function should return a single numeric value")
  }

  # 将其转换为字符串
  comparator_str <- rlang::as_string(rlang::get_expr(comparator.study))
  ipds_param_str <- rlang::as_string(rlang::get_expr(ipds.param))

  # Percentile/Basic/ CIs
  boot_ci_basic <- boot::boot.ci(boot.out = result_bootstraps,
                                 index = 1,
                                 type = "basic",
                                 conf = CIw)

  boot_ci_perc <- boot::boot.ci(boot.out = result_bootstraps,
                                index = 1,
                                type = "perc",
                                conf = CIw)

  boot_ci_bac <- boot::boot.ci(boot.out = result_bootstraps,
                               index = 1,
                               type = "bca",
                               conf = CIw)

  result_basic <- data.frame(
    comparator = comparator_str,
    param = ipds_param_str,
    type = "un-anchored",
    outcome = "bootstrap_Basic",
    effects = sprintf(paste0("%.", digits, "f"), boot_ci_basic$t0)
  )

  result_perc <- data.frame(
    comparator = comparator_str,
    param = ipds_param_str,
    type = "un-anchored",
    outcome = "bootstrap_Perc",
    effects = sprintf(paste0("%.", digits, "f"), boot_ci_perc$t0)
  )

  result_bac <- data.frame(
    comparator = comparator_str,
    param = ipds_param_str,
    type = "un-anchored",
    outcome = "bootstrap_BAC",
    effects = sprintf(paste0("%.", digits, "f"), boot_ci_bac$t0)
  )

  # 仅当 boot_ci 相应的列存在时添加列
  if (!is.null(boot_ci_basic$basic)) {
    result_basic[[lower_ci]] <- sprintf(paste0("%.", digits, "f"),
                                        boot_ci_basic$basic[4])
    result_basic[[upper_ci]] <- sprintf(paste0("%.", digits, "f"),
                                        boot_ci_basic$basic[5])
  } else {
    result_basic[[lower_ci]] <- NA
    result_basic[[upper_ci]] <- NA
  }

  if (!is.null(boot_ci_perc$percent)) {
    result_perc[[lower_ci]] <- sprintf(paste0("%.", digits, "f"),
                                       boot_ci_perc$percent[4])
    result_perc[[upper_ci]] <- sprintf(paste0("%.", digits, "f"),
                                       boot_ci_perc$percent[5])
  } else {
    result_perc[[lower_ci]] <- NA
    result_perc[[upper_ci]] <- NA
  }

  if (!is.null(boot_ci_bac$bca)) {
    result_bac[[lower_ci]] <- sprintf(paste0("%.", digits, "f"),
                                      boot_ci_bac$bca[4])
    result_bac[[upper_ci]] <- sprintf(paste0("%.", digits, "f"),
                                      boot_ci_bac$bca[5])
  } else {
    result_bac[[lower_ci]] <- NA
    result_bac[[upper_ci]] <- NA
  }

  results <- bind_rows(result_basic, result_perc, result_bac)

  # 绘图
  graphics::hist(result_bootstraps$t, main = "", xlab = "Boostrapped HR")
  graphics::abline(v = stats::quantile(result_bootstraps$t,
                                       probs = c(0.025, 0.5, 0.975)), lty = 2)

  plot <- grDevices::recordPlot()

  return(list(results = results, plot = plot))
}
