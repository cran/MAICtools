\name{AgD_eff}
\alias{AgD_eff}
\docType{data}
\title{Description of AgD_eff dataset}
\description{
  An example data frame containing aggregate results data from the comparator study.
}
\usage{
  data(AgD_eff)
}
\format{
  A data frame with X rows and Y variables:
  \describe{
    \item{STUDY}{Label of the comparator study, e.g., "Study XX-1".}
    \item{PARAM}{Subsets to be analyzed, e.g., "PFSINV", "OS".}
    \item{EST}{Point estimate of the effect size.}
    \item{CIL}{The lower confidence limit of the point estimate of the effect size.}
    \item{CIU}{The upper confidence limit of the point estimate of the effect size.}
  }
}
\examples{
  data(AgD_eff)
  head(AgD_eff)
}
