\name{pseudo}
\alias{pseudo}
\docType{data}
\title{Description of pseudo dataset}
\description{
  An example data frame containing pseudo patient data from the comparator study
}
\usage{
  data(pseudo)
}
\format{
  A data frame with X rows and Y variables:
  \describe{
    \item{SUBJID}{Subject Unique Identifier.}
    \item{PARAMCD}{Subsets to be analyzed, e.g., "PFSINV", "OS".}
    \item{ARM}{Label of the comparator study, = "Comparator".}
    \item{AVAL}{Survival or follow up time.}
    \item{CNSR}{The status indicator, 0 = event, 1 = censored.}
    \item{wt}{Weights, = 1.}
  }
}
\examples{
  data(pseudo)
  head(pseudo)
}
