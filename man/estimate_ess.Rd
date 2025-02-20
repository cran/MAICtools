% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_ess.R
\name{estimate_ess}
\alias{estimate_ess}
\title{Estimate Effective Sample Size (ESS)}
\arguments{
\item{ipds_wts}{A data frame containing individual patient data from the
intervention study, with a column containing the estimated weights (derived
using \code{\link{estimate_weights}}).}

\item{agds}{A data frame containing aggregate summary data from the
comparator study.}

\item{intervention.arm}{The name of the grouping column in the data frame
specified by \emph{ipds}, e.g., intervention.arm = TRT. The default is TRT.}

\item{comparator}{The name of the study column in the data frame specified
by \emph{agds}, e.g., comparator = STUDY. The default is STUDY.}

\item{comparator.study}{A character specifying the comparator study, which
must be quoted and cannot be empty (e.g., comparator.study = "Study XX-1").
This is the value of the study column in \emph{agds} set by the \emph{comparator}
parameter.}

\item{comparator.arm}{The name of the grouping column in the data frame
specified by \emph{agds}, e.g., comparator.arm = TRT. The default is TRT.}

\item{comparator.n}{A The name of the subjects number column in the data
frame specified by \emph{agds}, e.g., comparator.n = N. The default is N.}

\item{wt.col}{The name of the estimated weights column in the data frame
specified by \emph{ipds_wts}. The default is wt.}

\item{digits}{Specify the number of decimal places for the output results.}
}
\value{
A data frame containing effective sample size (ESS) after weighting.
}
\description{
Estimate Effective Sample Size (ESS)
}
\examples{
\donttest{
ess <- estimate_ess(
  ipds_wts = pts, agds = AgD_bl,
  intervention.arm = TRT,
  comparator = STUDY, comparator.study = "Study XX-1", comparator.arm = TRT,
  comparator.n = N)
ess
}

}
