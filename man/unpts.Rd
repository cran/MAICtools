\name{unpts}
\alias{unpts}
\docType{data}
\title{Description of unpts dataset}
\description{
  A combined data frame containing individual efficacy data from the intervention study and pseudo efficacy data from the comparator study.
}
\usage{
  data(unpts)
}
\format{
  A data frame with X rows and Y variables:
  \describe{
    \item{SUBJID}{Subject Unique Identifier.}
    \item{PARAMCD}{Subsets to be analyzed, e.g., "PFSINV", "OS".}
    \item{ARM}{Label of the study, "Intervention" for the intervention study and "Cmparator" for the comparator study.}
    \item{AVAL}{Survival or follow up time.}
    \item{CNSR}{The status indicator, 0 = event, 1 = censored.}
    \item{wt}{Weights to be used.}
  }
}
\examples{
  data(unpts)
  head(unpts)
}
