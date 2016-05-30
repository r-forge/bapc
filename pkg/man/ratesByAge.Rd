\name{ratesByAge}
\alias{ratesByAge}
\title{
    Function to plot age-specific incidence or mortality rates along time
}
\description{
    This function is based on the rateplot function of the Epi package. 
    It uses a APCList object and plots age-specific rates 
    against period.
}
\usage{
ratesByAge(APCList, age, per, scale=10^3, log.ax = "", ann = TRUE, grid = TRUE, ...)
}
\arguments{
  \item{APCList}{
    A \code{APCList} object
}
  \item{scale}{
    Defining the scale of the rates.
}
  \item{age}{
    Numerical vector giving the means of the age-classes.
}
  \item{per}{
    Numerical vector giving the means of the periods. 
}
  \item{log.ax}{
    Character where '"y"' and '"r"' both refer to the rate scale. Defaults to '""'.
}
  \item{ann}{
    Should the curves be annotated? Default: TRUE.
}
  \item{grid}{
    Logical indicating whether a background grid should be drawn. Default: TRUE.
}
  \item{\dots}{
    Additional arguments pssed on to 'rateplot'.
}
}
\details{

     Zero values of the rates are ignored. They are not shown in the plot.
}
\value{
     'NULL'. The function is used for its side-effect, the plot.
}
\references{
  Bendix Carstensen, Martyn Plummer, Esa Laara, Michael Hills (2012).
  Epi: A Package for Statistical Analysis in Epidemiology. R package
  version 1.1.40. URL http://CRAN.R-project.org/package=Ep
}
\author{
    Andrea Riebler
}

\seealso{
    \code{\link{rateplot}}
}
\examples{
   data(FemLCSweden)
   data(FemPYSweden)

   lc_sweden <- APCList(FemLCSweden, FemPYSweden, gf=5)

   ratesByAge(lc_sweden, 10^3, age=seq(27,82,by=5), per=1951:2010, col=topo.colors(nage(lc_sweden)))
}
