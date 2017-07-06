\name{BAPC}
\alias{BAPC}
\title{
    Function to project age-specific mortality or disease rates using INLA
}
\description{
    A Bayesian age-period-cohort model fitted using a Poisson model within INLA is used to project
    mortality or disease rates. Age, period and/or cohort effects are either modelled using a random walk
    of second order (RW2), or fixed effect (drift).
}
\usage{
  BAPC(APCList, predict = list(npredict = 0, retro = TRUE),
  model = list(age = list(model = "rw2", prior="loggamma", param = c(1, 0.00005), initial = 4, scale.model=FALSE),
    period = list(include = TRUE, model = "rw2", prior="loggamma", param = c(1, 0.00005), initial = 4, scale.model=FALSE),
    cohort = list(include = TRUE, model = "rw2", prior="loggamma", param = c(1, 0.00005), initial = 4, scale.model=FALSE),
    overdis = list(include = TRUE, model = "iid", prior="loggamma", param = c(1, 0.005), initial = 4)),
    secondDiff = FALSE, stdweight = NULL, verbose = FALSE)
}
\arguments{
  \item{APCList}{
    a \code{APCList} object
}
  \item{predict}{
    a list specifying how many periods are to be projected and whether existing counts shoud be removed and projected.
    The first argument \code{npredict} sets the number of periods to predict. The second argument \code{retro} is a boolean
    indicating whether the past \code{npredict} periods should be removed or predicted or whether future rates are to be predicted.
    If \code{retro} is FALSE data cells for which the corresponing observation count is set to NA will be projected.}
  \item{model}{
    a list composed of four arguments: age, period, cohort, overdis. For each argument a separate list is to be specified defining whether the component should be included and if so, which model should be used. Possible arguments are:
    \describe{
        \item{include}{A Boolean indictating whether the component should be included in the model. Age is per default always included.}
        \item{model}{A character indicating the model to be used. Possible choices are \code{rw2} and \code{drift}.}
        \item{prior}{Prior distribution for the (log) precision parameter, such as "loggamma" or "pc.prec" for the penalised-complexity prior.}
        \item{param}{Parameters for the prior distribution for the precision parameter.}
        \item{initial}{Initial value for the precision parameter (on log-scale). }
        \item{scale.model}{Logical. If TRUE then scale the RW1 and RW2 models so that their (generlized) variance is 1, see Sorbye and Rue (2014).}}
    }
  \item{secondDiff}{
    Boolean (default:FALSE) indicating whether summary estimates for the second differences of age, period and cohort effects should be computed.
  }
  \item{stdweight}{
    Numeric vector with length equal to the number of age groups used to derive age-standardized projections.
    If the weights do not sum to one, they will be normalised internally. If no weights are provided, there will 
    no age-standardized projections be computed.
    }
  \item{verbose}{
    Boolean (default:FALSE) indicating whether the program should run in a verbose model.
    }    
}
\value{
    An \code{APCList} object.
}
\references{
  Havard Rue, Sara Martino, and Nicholas Chopin (2009). Approximate
  Bayesian Inference for Latent Gaussian Models Using Integrated Nested
  Laplace Approximations. Journal of the Royal Statistical Society B,
  71, 319-392. (www.r-inla.org)


  Sigrunn Holbek Sorbye and Havard Rue (2014). "Scaling intrinsic Gaussian Markov random field priors in spatial modelling." Spatial Statistics 8: 39-51.
}
\author{
    Andrea Riebler and Leonhard Held
}

\seealso{
    \code{inla}
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
      secondDiff=FALSE, stdweight=whostandard[6:17,2], verbose=FALSE)
}
}
}
