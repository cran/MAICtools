#' @title Functions for the Estimation of Propensity Weights
#'
#' @param ipds A data frame containing individual patient data from the
#' intervention study, with baseline characteristic variables for matching.
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
#' @param opt.method The optim method to be used. The default is "BFGS".
#' @param seed The seed for centralized variable missing value imputation
#' (KNN method).
#' @param ... Refer to \link[stats:optim]{optim} for additional parameters.
#'
#' @importFrom rlang enquo as_name
#' @importFrom assertthat assert_that
#' @importFrom dplyr union filter mutate select ends_with
#' @importFrom magrittr %>%
#' @importFrom stats optim
#' @importFrom VIM kNN
#' @importFrom purrr map_dfr
#' @importFrom data.table :=
#'
#' @return A data frame containing individual patient data, calculated weights,
#' and rescaled weights.
#' @export
#'
#' @examples
#' \donttest{
#' cov <- list(
#'   c("ECOG", "SMK", "METBRAIN"),
#'   c("BMI", "DIAG")
#' )
#'
#' pts <- estimate_weights(
#'   ipds = IPD,
#'   agds = AgD_bl,
#'   matching.list = cov,
#'   intervention.arm = TRT,
#'   comparator = STUDY,
#'   comparator.study = "Study XX-1",
#'   comparator.arm = TRT
#' )
#' }
#' @name estimate_weights
NULL

# 声明全局变量
utils::globalVariables(c("TRT", "STUDY", "wt", ".data"))

estimate_weights <- function(ipds, agds, matching.list,
                             intervention.arm = TRT,
                             comparator = STUDY, comparator.study,
                             comparator.arm = TRT,
                             opt.method = "BFGS", seed = 123456, ...) {

  # 使用enquo捕获参数符号
  intervention.arm <- rlang::enquo(intervention.arm)
  comparator <- rlang::enquo(comparator)
  comparator.arm <- rlang::enquo(comparator.arm)

  # checking--------------------------------------------------------------------
  # 检查ipds和agds是否为data.frame类型
  assertthat::assert_that(
    is.data.frame(ipds),
    msg = "'ipds' is expected to be a data frame")
  assertthat::assert_that(
    is.data.frame(agds),
    msg = "'agds' is expected to be a data frame")

  # 检查intervention.arm和comparator.arm是否在ipds和agds中
  assertthat::assert_that(
    rlang::as_name(intervention.arm) %in% names(ipds),
    msg = paste0(rlang::as_name(intervention.arm),
                 " can not be found in ipds"))
  assertthat::assert_that(
    rlang::as_name(comparator.arm) %in% names(agds),
    msg = paste0(rlang::as_name(comparator.arm),
                 " can not be found in agds"))

  # 检查matching.list的类型及其变量是否在数据框中
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

  # 查找连续变量的mean和sd是否在 agds 中
  if (length(matching.list[[2]]) > 0) {
    match_cov_cont_mean <- paste0(unlist(matching.list[[2]]), ".mean")
    missing_vars_cont_mean <-
      match_cov_cont_mean[!match_cov_cont_mean %in% names(agds)]
    if (length(missing_vars_cont_mean) > 0) {
      stop(paste(missing_vars_cont_mean, collapse = "; "),
           " cannot be found in agds")
    }

    match_cov_cont_sd <- paste0(unlist(matching.list[[2]]), ".sd")
    missing_vars_cont_sd <-
      match_cov_cont_sd[!match_cov_cont_sd %in% names(agds)]
    if (length(missing_vars_cont_sd) > 0) {
      stop(paste(missing_vars_cont_sd, collapse = "; "),
           " cannot be found in agds")
    }
  }

  # 定义必要的函数
  # Objective function
  objfn <- function(a1, X) {sum(exp(X %*% a1))}

  # Gradient function
  gradfn <- function(a1, X) {colSums(sweep(X, 1, exp(X %*% a1), "*"))}

  # 保存intervention.arm类型和顺序 ---------------------------------------------
  intervention_levels <- levels(ipds[[rlang::as_name(intervention.arm)]])
  intervention_arm_as_fac <- is.factor(ipds[[rlang::as_name(intervention.arm)]])

  # matching--------------------------------------------------------------------
  match_cov_centered <- c(
    unlist(lapply(matching.list[[1]], function(x) paste0(x, "_centered"))),
    unlist(lapply(matching.list[[2]], function(x) paste0(x, ".mean_centered"))),
    unlist(lapply(matching.list[[2]], function(x) paste0(x, ".sd_centered")))
  )

  intervention_values <- unique(ipds[[rlang::as_name(intervention.arm)]])

  # 定义函数来处理中心化和权重
  process_weights <- function(intervention_value) {

    intervention_data <- ipds %>%
      dplyr::filter(!!intervention.arm == intervention_value)
    intervention_data_centered <- intervention_data
    comparator_data <- agds %>%
      dplyr::filter(!!comparator == comparator.study
                    & !!comparator.arm == intervention_value)

    # 如果类别指标的向量不为空，则进行二分类变量中心化处理
    if (length(matching.list[[1]]) > 0) {
      for (i in matching.list[[1]]) {
        value <- comparator_data[[i]]

        if (is.numeric(value) && length(value) == 1) {
          intervention_data_centered <- intervention_data_centered %>%
            dplyr::mutate(!!paste0(i, "_centered") := .data[[i]] - value)
        } else {
          stop(paste("variable ", i,
                     " in agds does not have a unique numeric value"))
        }
      }
    }

    # 如果连续指标的向量不为空，则进行连续变量中心化处理
    if (length(matching.list[[2]]) > 0) {
      for (i in matching.list[[2]]) {
        mean_val <- comparator_data[[paste0(i, ".mean")]]
        sd_val <- comparator_data[[paste0(i, ".sd")]]

        if (is.numeric(mean_val)
            && length(mean_val) == 1
            && is.numeric(sd_val)
            && length(sd_val) == 1) {
          intervention_data_centered <- intervention_data_centered %>%
            dplyr::mutate(
              !!paste0(i, ".mean_centered") := .data[[i]] - mean_val,
              !!paste0(i, ".sd_centered")
              := (.data[[i]]^2) - (mean_val^2 + sd_val^2))
        } else {
          stop("variable ", i,
               " in agds does not have unique numeric mean or sd.")
        }
      }
    }

    # 检查中心化变量缺失比例并用KNN处理
    columun_with_centered <- grep("_centered$",
                                  colnames(intervention_data_centered),
                                  value = TRUE)

    set.seed(seed)
    for (i in columun_with_centered) {

      missing_ratio <- mean(is.na(intervention_data_centered[[i]]))

      if (missing_ratio > 0 && missing_ratio <= 0.15) {
        intervention_data_centered[[i]] <- VIM::kNN(
          intervention_data_centered,
          variable = i, k = 5)[[i]]
      }
      else if (missing_ratio > 0.15) {
        warning(paste0("Too many missing values for ", i))
      }
    }

    # Optimise Q(b) using Newton-Raphson techniques
    opt1 <- stats::optim(
      par = rep(0, dim(
        as.data.frame(intervention_data_centered[, match_cov_centered]))[2]),
      fn = objfn,
      gr = gradfn,
      X = as.matrix(intervention_data_centered[, match_cov_centered]),
      method = opt.method,
      ...
    )

    # Calculate weights for intervention data and combine with dataset
    wts <- dplyr::mutate(
      intervention_data_centered,
      wt = as.vector(exp(as.matrix(
        intervention_data_centered[, match_cov_centered]) %*% opt1$par)),
      # rescaled weights
      wt_rs = (wt / sum(wt)) * nrow(intervention_data_centered),
    )

    return(wts)
  }

  data_with_wts <- purrr::map_dfr(intervention_values, process_weights) %>%
    dplyr::select(-dplyr::ends_with("_centered"))

  # 恢复intervention.arm的原始类型和顺序
  if (intervention_arm_as_fac) {
    data_with_wts <- data_with_wts %>%
      mutate(!!rlang::as_name(intervention.arm)
             := factor(.data[[rlang::as_name(intervention.arm)]],
                       levels = intervention_levels)) %>%
      arrange(!!intervention.arm)
  }

  return(data_with_wts)
}
