\name{pts}
\alias{pts}
\docType{data}
\title{Description of pts dataset}
\description{
  An example data frame containing individual patient data and estimated weights.
}
\usage{
  data(pts)
}
\format{
  A data frame with X rows and Y variables:
  \describe{
    \item{SUBJID}{Subject Unique Identifier.}
    \item{TRT}{Grouping variable, e.g., "active" or "control".}
    \item{BPDL1}{Stratification factors for stratified analysis.}
    \item{CNSBRAIN}{Stratification factors for stratified analysis.}
    \item{AGEGR}{Stratification factors for stratified analysis.}
    \item{WEIGHT}{Baseline characteristic variables for matching or summarizing.}
    \item{HEIGHT}{Baseline characteristic variables for matching or summarizing.}
    \item{BMI}{Baseline characteristic variables for matching or summarizing.}
    \item{DIAG}{Baseline characteristic variables for matching or summarizing.}
    \item{SEX}{Baseline characteristic variables for matching or summarizing.}
    \item{ECOG}{Baseline characteristic variables for matching or summarizing.}
    \item{SMK}{Baseline characteristics variables for matching.}
    \item{METBRAIN}{Baseline characteristics variables for matching.}
    \item{METLIVER}{Baseline characteristics variables for matching.}
    \item{PARAMCD}{Subsets to be analyzed, e.g., "PFSINV", "OS".}
    \item{AVAL}{Survival or follow up time.}
    \item{CNSR}{The status indicator, 0 = event, 1 = censored.}
    \item{RESP}{Response status, 1 = responder, 0 = non-responder.}
    \item{wt}{Estimated propensity weights.}
    \item{wt_rs}{Estimated rescaled propensity weights.}
  }
}
\examples{
  data(pts)
  head(pts)
}
