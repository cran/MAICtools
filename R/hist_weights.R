#' @title Histograms of Weights and Rescaled Weights Distributions
#'
#' @param ipds_wts A data frame containing individual patient data from the
#' intervention study, with a column containing the estimated weights
#' (derived using \code{\link{estimate_weights}}).
#' @param intervention.arm The name of the grouping column in the data frame
#' specified by *ipds*, e.g., `intervention.arm = TRT`. The default is `TRT`.
#' @param wt.col The name of the estimated weights column in the data frame
#' specified by *ipds_wts*. The default is `wt`.
#' @param rswt.col The name of the estimated rescaled weights column in the
#' data frame specified by *ipds_wts*. The default is `wt_rs`.
#' @param bin The number of bins or bars of the histogram.
#' @param xstepby An integer guiding the breaks on the X-axis.
#' @param ystepby An integer guiding the breaks on the Y-axis.
#' @param ... Refer to \link[ggplot2:geom_histogram]{geom_histogram} for
#' additional parameters.
#'
#' @importFrom rlang enquo as_name
#' @importFrom assertthat assert_that
#' @importFrom dplyr select rename
#' @importFrom magrittr %>%
#' @importFrom graphics hist
#' @importFrom ggplot2 ggplot aes geom_histogram facet_grid scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous scale_fill_manual theme_bw theme
#' @importFrom ggplot2 element_text element_line element_blank unit labs
#' @importFrom tidyr pivot_wider
#'
#' @return Histograms of weights and rescaled weights distributions.
#' @export
#'
#' @examples
#' \donttest{
#' hist_weights(pts, intervention.arm = TRT, xstepby = 2, ystepby = 50)
#' }
#' @name hist_weights

# 声明全局变量
utils::globalVariables(c("TRT", "wt", "wt_rs", "value"))

hist_weights <- function(ipds_wts,
                         intervention.arm = TRT,
                         wt.col = wt,
                         rswt.col = wt_rs,
                         bin = 30, xstepby = 1, ystepby = 30, ...) {

  # 使用enquo捕获参数符号
  intervention.arm <- rlang::enquo(intervention.arm)
  wt.col <- rlang::enquo(wt.col)
  rswt.col <- rlang::enquo(rswt.col)

  # checking--------------------------------------------------------------------
  # 检查ipds是否为data.frame类型
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

  # histogram
  wt_data <- ipds_wts %>%
    dplyr::select(!!intervention.arm, !!wt.col, !!rswt.col) %>%
    dplyr::rename(Intervention = !!intervention.arm,
                  Weights = !!wt.col, `Rescaled Weights` = !!rswt.col) %>%
    tidyr::pivot_longer(cols = c("Weights", "Rescaled Weights"),
                        names_to = "key", values_to = "value")

  wt_data$key <- factor(wt_data$key, levels = c("Weights", "Rescaled Weights"))

  max_count <- max(vapply(levels(wt_data$key), function(k) {
    max(graphics::hist(wt_data$value[wt_data$key == k], plot = FALSE,
                       breaks = bin)$counts)
  }, numeric(1)))

  hist_plot <- ggplot2::ggplot(wt_data,
                               ggplot2::aes(x = value, fill = "#595959")) +
    ggplot2::geom_histogram(bins = bin, position = "dodge", ...) +
    ggplot2::facet_grid(Intervention ~ key) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, max(wt_data$value, na.rm = TRUE), by = xstepby)) +
    ggplot2::scale_y_continuous(
      breaks = seq(0, max_count, by = ystepby)) +
    ggplot2::scale_fill_manual(values = "#595959") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 16),
      axis.text = ggplot2::element_text(size = 12),
      strip.text.x = ggplot2::element_text(size = 12),
      strip.text.y = ggplot2::element_text(size = 14),
      legend.position = "none",
      panel.grid.major = ggplot2::element_line(linewidth = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      panel.spacing = grid::unit(1, "lines")
    ) +
    ggplot2::labs(y = "Frequency", x = "Weight")

  return(hist_plot)
}
