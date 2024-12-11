\name{AgD_bl}
\alias{AgD_bl}
\docType{data}
\title{Description of AgD_bl dataset}
\description{
  An example data frame containing aggregate summary data from the comparator study.
}
\usage{
  data(AgD_bl)
}
\format{
  A data frame with X rows and Y variables:
  \describe{
    \item{STUDY}{Label of the comparator study, e.g., "Study XX-1".}
    \item{TRT}{Grouping variable, e.g., "active" or "control".}
    \item{N}{Number of subjects in each group.}
    \item{AGEGR}{Baseline characteristics variables for matching.}
    \item{SEX}{Baseline characteristics variables for matching.}
    \item{ECOG}{Baseline characteristics variables for matching.}
    \item{SMK}{Baseline characteristics variables for matching.}
    \item{METBRAIN}{Baseline characteristics variables for matching.}
    \item{METLIVER}{Baseline characteristics variables for matching.}
    \item{BMI.mean}{Baseline characteristics variables for matching.}
    \item{BMI.sd}{Baseline characteristics variables for matching.}
    \item{DIAG.mean}{Baseline characteristics variables for matching.}
    \item{DIAG.sd}{Baseline characteristics variables for matching.}
  }
}
\examples{
  data(AgD_bl)
  head(AgD_bl)
}
