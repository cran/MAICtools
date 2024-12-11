#' @title Generate a Kaplan-Meier Plot with Individual Efficacy Data and Pseudo
#' Efficacy Data.
#'
#' @param unds_wts A combined data frame containing individual efficacy data
#' from the intervention study and pseudo efficacy data from the comparator
#' study.
#' @param unds.arm The name of the grouping column in the combined data frame
#' specified by *unds_wts*, e.g., comparator.arm = TRT. The default is TRT.
#' @param unds.param.var The name of the column that specifies only a subset
#' of the rows of the data to be used.
#' @param unds.param A character specifying the subset of the rows to be used.
#' This is the value of the column set by the *unds.param.var*.
#' @param time The name of the survival or follow up time column in the
#' combined data frame.
#' @param status The status indicator, normally 0 = event, 1 = censored. Can
#' be reseted using the *event* parameter.
#' @param event A numeric value that represents the survival status, 0 = event,
#' 1 = censored.
#' @param wt.col The name of the estimated weights column in the data frame
#' specified by *unds_wts*. The default is wt.
#' @param km.xlim A numeric value specifying the right limit of the scale on
#' the X-axis.
#' @param xstepby An integer guiding the breaks on the X-axis.
#' @param km.ylim A numeric value specifying the upper limit of the scale on
#' the Y-axis.
#' @param ystepby An integer guiding the breaks on the Y-axis.
#' @param xlab A character giving label of the X-axis. The default is
#' "Time (Months)".
#' @param ylab A character giving label of the Y-axis. The default is
#' "Survival probability".
#' @param km.legend A character vector of length >=1 to appear in the legend.
#' @param km.title A character used to set the main title at the top.
#' @param ... Refer to \link[survminer:ggsurvplot]{ggsurvplot} for additional
#' parameters..
#'
#' @importFrom rlang enquo as_name get_expr
#' @importFrom assertthat assert_that
#' @importFrom dplyr filter
#' @importFrom survival Surv survfit
#' @importFrom survminer ggsurvplot theme_cleantable
#' @importFrom ggplot2 theme margin element_text element_blank
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous annotate
#'
#' @return A Kaplan-Meier plot object that contains individual efficacy data
#' from the intervention study and pseudo efficacy data from the comparator
#' study.
#' @export
#'
#' @examples
#' \donttest{
#' unanchored_kmplot(
#'   unds_wts = unpts, unds.arm = ARM,
#'   unds.param.var = PARAMCD, unds.param = "OS",
#'   time = AVAL, status = CNSR, event = 0,
#'   wt.col = wt, km.xlim = 35, xstepby = 3,
#'   km.legend = c("Arm A", "ARM B"),
#'   km.title = "AAAA")
#' }
#'
#' @name unanchored_kmplot
NULL

# 声明全局变量
utils::globalVariables(c("ARM", "wt"))

unanchored_kmplot <- function(unds_wts, unds.arm = ARM,
                              unds.param.var = NULL, unds.param = NULL,
                              time = NULL, status = NULL, event = 0,
                              wt.col = wt,
                              km.xlim, xstepby, km.ylim = 1, ystepby = 0.25,
                              xlab = "Time (Months)",
                              ylab = "Survival probability",
                              km.legend = NULL, km.title = NULL, ...) {
  # 使用 enquo 捕获参数符号
  unds.arm <- rlang::enquo(unds.arm)
  unds.param.var <- rlang::enquo(unds.param.var)
  time <- rlang::enquo(time)
  status <- rlang::enquo(status)
  wt.col <- rlang::enquo(wt.col)

  # checking--------------------------------------------------------------------
  # 检查unds_wts是否为data.frame类型并删除!<-参考
  assertthat::assert_that(
    is.data.frame(unds_wts),
    msg = "'unds_wts' is expected to be a data frame")
  # 检查变量unds.arm、wt.col、time、status是否在数据集unds_wts中
  assertthat::assert_that(
    rlang::as_name(unds.arm) %in% names(unds_wts),
    msg = paste0(rlang::as_name(unds.arm),
                 " can not be found in unds_wts"))
  assertthat::assert_that(
    rlang::as_name(wt.col) %in% names(unds_wts),
    msg = paste0(rlang::as_name(wt.col),
                 " can not be found in unds_wts"))
  assertthat::assert_that(
    rlang::as_name(time) %in% names(unds_wts),
    msg = paste0(rlang::as_name(time),
                 " can not be found unds_wts"))
  assertthat::assert_that(
    rlang::as_name(status) %in% names(unds_wts),
    msg = paste0(rlang::as_name(status),
                 " can not be found in unds_wts"))

  # 处理参数变量筛选
  if (!is.null(unds.param) && !is.null(rlang::get_expr(unds.param.var))) {
    km_ds <- dplyr::filter(!!unds.param.var == unds.param, .data = unds_wts)
  } else {
    km_ds <- unds_wts
  }

  #处理legend_labels
  legend_labels <- if (!is.null(km.legend)) {
    km.legend
  } else {
    levels(unds_wts$ARM)
  }

  # 计算 x轴最大值并处理y轴符号边界
  handle_xlim <- function(unds_wts) {
    max_val <- ceiling(max(unds_wts[[rlang::as_name(time)]], na.rm = TRUE))
    time_step <- 0

    if (max_val <= 24) {
      time_step <- 2
    } else {
      time_step <- 3
    }

    return(list(time_max = max_val, time_step = time_step))
  }

  params <- handle_xlim(unds_wts)
  km.xlim_real <- if (is.null(km.xlim)) params$time_max else km.xlim
  xstepby_real <- if (is.null(xstepby)) params$time_step else xstepby

  # 构建公式
  formula <- stats::as.formula(
    paste0("survival::Surv(", rlang::as_name(time), ", ",
           rlang::as_name(status), " == ", event, ") ~ ",
           rlang::as_name(unds.arm)))

  # 构建 KM 模型
  km_list <- do.call(survival::survfit,
                     list(formula = formula,
                          data = unds_wts,
                          weights = rlang::eval_tidy(wt.col, unds_wts),
                          type = "kaplan-meier", conf.type = "log-log"))

  # 绘制 KM 曲线
  km_plot <- survminer::ggsurvplot(
    km_list, data = km_ds, surv.median.line = "hv",
    risk.table = FALSE, tables.theme = survminer::theme_cleantable(),
    palette = c("#CD5C5C", "#00468B"), xlab = xlab, ylab = ylab,
    xlim = c(0, km.xlim_real), break.time.by = xstepby_real,
    font.x = 14, legend = c(0.8, 0.8), legend.labs = legend_labels,
    title = km.title, font.title = c(14), ...)

  median <- paste0(sprintf(paste0("%.", 1, "f"), summary(km_list)$table[13]),
                   " months vs. ",
                   sprintf(paste0("%.", 1, "f"), summary(km_list)$table[14]),
                   " months")

  # 添加绘图
  km_plot$plot <- km_plot$plot + ggplot2::scale_x_continuous(
    expand = c(0, 0), limits = c(0, km.xlim_real),
    breaks = seq(0, km.xlim_real, xstepby_real)) +
    ggplot2::scale_y_continuous(
      expand = c(0, 0),
      limits = c(0, km.ylim),
      breaks = seq(0, km.ylim, ystepby)) +
    ggplot2::annotate(
      geom = "text",
      x = km.xlim_real / 6.5, y = 0.53,
      label = median,
      size = 4.5) +
    ggplot2::theme(
      plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10),
      axis.text.x = ggplot2::element_text(size = 14),
      axis.text.y = ggplot2::element_text(size = 14),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 14),
      legend.direction = "vertical"
    )

  return(km_plot$plot)
}
