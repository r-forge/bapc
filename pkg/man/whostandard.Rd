\name{whostandard}
\alias{whostandard}
\docType{data}
\title{
    WHO world standard population
}
\description{
    Percentage of the population in each 5-year age group in the WHO World Standard population from zero to 85+ years.
}
\usage{data(whostandard)}
\format{
  A data frame with 18 observations on the following 2 variables.
  \describe{
    \item{\code{agegroup}}{a factor with levels \code{0-4} \code{10-14} \code{15-19} \code{20-24} \code{25-29} \code{30-34} \code{35-39} \code{40-44} \code{45-49} \code{50-54} \code{55-59} \code{5-9} \code{60-64} \code{65-69} \code{70-74} \code{75-79} \code{80-84} \code{85+}}
    \item{\code{whoStandard}}{a numeric vector containing the percentage of the population in the different age groups.}
  }
}
\source{
Ahmad, O., Boschi-Pinto, C., Lopez, A. D., Murray, C. J. L., Lozano, R. and Inoue, M.
(2001). Age standardization of rates: a new WHO standard, GPE Discussion Paper
Series, Geneva. (GPE Discussion Paper No. 31).
}
\examples{
data(whostandard)
}
\keyword{datasets}
