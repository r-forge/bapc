\name{summaryHyper}
\alias{summaryHyper}
\title{
    Function to get summary information for the hyperparamters, i.e. precision or variance parameters, in the model.
    }
\description{
    This function returns summary information for the hyperparameters. The default is to use precision parameters, but summary information for variance parameters is available as well. 
}
\usage{
  summaryHyper(APCList, var=TRUE)
}
\arguments{
  \item{APCList}{
    a \code{APCList} object
   }
  \item{var}{
    Boolean (default:TRUE) indicating whether the summary information should be returned on the variance parameter scale of precision parameter scale.
}
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
  
   summaryHyper(result, var=FALSE)
}
}
}
