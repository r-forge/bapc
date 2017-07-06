\name{summarySecDiff}
\alias{summarySecDiff}
\title{
    Function to get summary information for the second differences of age, period or cohort effects.
}
\description{
    Due to the identifiability problem of the APC model are the latent effects not interpretable, however, the second differences representing measure of curvature are. This function returns summary information for the second differences either on the internal log scale or the exponential scale for age, period or cohort effects. On the exponential scale, these contrasts represent ratios of two adjacent relative
risks. Caution, in order to extract the summary information the option \code{secondDiff} needs to be set in the function \code{BAPC}.
}
\usage{
  summarySecDiff(APCList, variable="age", log=FALSE)
}
\arguments{
  \item{APCList}{
    a \code{APCList} object
  }
  \item{variable}{
    a string indicating for which of "age", "period" or "cohort" effects, the summary information is desired.
  }
  \item{log}{
    Boolean (default:FALSE) indicating whether the summary information should be returned on the log scale, i.e. the   real second differences of age effects, say, or for the exponential transformation to get in on the observation scale.}
    }
\value{
    A matrix including mean, standard deviation and quantile information.
}
\author{
    Andrea Riebler
}

\seealso{
    \code{BAPC}
}
\examples{
\dontrun{
if(requireNamespace("INLA", quietly = TRUE)) {
   require(INLA)
   data(FemLCSweden)
   data(FemPYSweden)
   data(whostandard)

   lc_sweden <- APCList(FemLCSweden, FemPYSweden, gf=5)

   result <- BAPC(lc_sweden, predict=list(npredict=10, retro=TRUE),
      secondDiff=TRUE, stdweight=NULL, verbose=FALSE)
  
   summarySecDiff(result, variable="period", log=FALSE)
}
}
}
