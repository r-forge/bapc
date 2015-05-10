\name{APCList}
\docType{class}
\alias{APCList}
\alias{APCList-class}
\alias{APCList,data.frame,data.frame,numeric-method}
\alias{[,APCList,ANY,ANY,ANY-method}
\alias{nage}
\alias{nage,APCList-method}
\alias{nperiod}
\alias{nperiod,APCList-method}
\alias{ncohort}
\alias{ncohort,APCList-method}
\alias{gridfactor}
\alias{gridfactor,APCList-method}
\alias{epi}
\alias{epi,APCList-method}
\alias{pyrs}
\alias{pyrs,APCList-method}
\alias{agelabels}
\alias{agelabels,APCList-method}
\alias{agelabels<-}
\alias{agelabels<-,APCList-method}
\alias{periodlabels}
\alias{periodlabels,APCList-method}
\alias{periodlabels<-}
\alias{periodlabels<-,APCList-method}
\alias{cohortlabels}
\alias{cohortlabels,APCList-method}
\alias{cohortlabels<-}
\alias{cohortlabels<-,APCList-method}
\alias{ageindex,APCList-method}
\alias{ageindex}
\alias{periodindex,APCList-method}
\alias{periodindex}
\alias{periodindex<-}
\alias{periodindex<-,APCList-method}
\alias{cohortindex,APCList-method}
\alias{cohortindex}
\alias{cohortindex<-}
\alias{cohortindex<-,APCList-method}
\alias{overdisindex,APCList-method}
\alias{overdisindex}
\alias{overdisindex<-}
\alias{overdisindex<-,APCList-method}
\alias{agespec.rate}
\alias{agespec.rate,APCList-method}
\alias{agespec.rate<-}
\alias{agespec.rate<-,APCList-method}
\alias{agespec.proj}
\alias{agespec.proj,APCList-method}
\alias{agespec.proj<-}
\alias{agespec.proj<-,APCList-method}
\alias{agestd.rate}
\alias{agestd.rate,APCList-method}
\alias{agestd.rate<-}
\alias{agestd.rate<-,APCList-method}
\alias{agestd.proj}
\alias{agestd.proj,APCList-method}
\alias{agestd.proj<-}
\alias{agestd.proj<-,APCList-method}
\alias{stdweight}
\alias{stdweight,APCList-method}
\alias{stdweight<-}
\alias{stdweight<-,APCList-method}
\alias{stdobs}
\alias{stdobs,APCList-method}
\alias{stdobs<-}
\alias{stdobs<-,APCList-method}
\alias{inlares}
\alias{inlares,APCList-method}
\alias{inlares<-}
\alias{inlares<-,APCList-method}

\title{Class \code{"APCList"}}
\description{
    This S4 class captures incidence or mortality data togehter with the 
    corresponding person-years and a gridfactor, indicating how much wider
    the agegroup intervals are compared to the period interval length.
}
\section{Constructor}{
    Creates a APCList object:
    \describe{
    \item{}{
        \code{APCList(epi, pyrs, gf, agelab=character(), periodlab=character(), cohortlab=character())}
        \describe{
            \item{\code{epi}}{A \code{data.frame} object of disease or mortality counts. 
                The number of rows are equal to the number of periods and the number of columns
                equal to the number of age groups. Data cells to be projected should be set to NA.}
            \item{\code{pyrs}}{A \code{data.frame} object of person years of population
                estimates. This data frame needs to have the same dimension as \code{epi}.
                Thus, the number of rows are equal to the number of periods and the number of columns
                equal to the number of age groups. Caution: Person years estimates for data cells that 
                are to be projected must still be provided.}
            \item{\code{gf}}{A grid factor indicating how much wider the agegroup intervals are compared
                to the period interval length. For example if periods are given annualy but age-groups
                are given for 5-year groups, the grid factor is equal to five.}
            \item{\code{agelab}}{A character vector with age group labels.  }
            \item{\code{periodlab}}{A character vector with period labels.}
            \item{\code{cohortlab}}{A character vector with cohort labels.}        
            \item{\code{ageindex}}{An integer index vector for age groups (derived internally).}  
            \item{\code{periodindex}}{An integer index vector for periods (derived internally).}
            \item{\code{cohortindex}}{An integer index vector for cohort (derived internally).}    
            \item{\code{overdisindex}}{An integer index vector for the overdispersion (derived internally).}    
            \item{\code{stdweight}}{Vector of length I giving one weight for each age group to compute age-standardized projections.}
            \item{\code{stdobs}}{Vector of length I containing age-standardized observations.}
            \item{\code{npred}}{Number of periods that were predicted.}
            \item{\code{agespec.rate}}{A list of length I, one for each age group. Each list element contains a matrix of dimension J * 4, where J denotes the number of periods. The columns present the 2.5\% CI, the mean of the rate, the 97.5\% CI and the standard deviation of the rate.}
            \item{\code{agespec.proj}}{A list of length I, one for each age group. Each list element contains a matrix of dimension J * 4, where J denotes the number of periods. The columns present the 2.5\% CI, the mean of the predictive distribution, the 97.5\% CI and the standard deviation of the predictive distribution.}
            \item{\code{agestd.rate}}{A matrix of dimension J * 4 containing 2.5\% CI, mean and the 97.5\% CI and the standard deviation of age standardized rates for all J periods.}
            \item{\code{agestd.proj}}{A matrix of dimension J * 4 containing 2.5\% CI, mean and the 97.5\% CI and the standard deviation of age standardized projections for all J periods.}
            \item{\code{inlares}}{A list containing the INLA results object.}}
    }
    }
}
\section{Methods}{
  \describe{
    \item{x[i,j]}{\code{signature(x = "APCList")}: Creates a APCList 
        object, keeping only the \code{i} periods and \code{j} age groups. }
    \item{nage}{\code{signature(object = "APCList")}: 
        get the number of agegroups.}
    \item{nperiod}{\code{signature(object = "APCList")}: 
        get the number of periods.}
    \item{ncohort}{\code{signature(object = "APCList")}: 
        get the number of cohorts.}
    \item{gridfactor}{\code{signature(object = "APCList")}: 
        get the gridfactor.}
    \item{epi}{\code{signature(object = "APCList")}: 
        get the disease/mortality counts.}
    \item{pyrs}{\code{signature(object = "APCList")}: 
        get the person years.}
    \item{ageindex}{\code{signature(x = "APCList")}: 
        get the \code{ageindex} used in INLA formula}
    \item{periodindex}{\code{signature(x = "APCList")}: 
        get the \code{periodindex} used in the INLA formula}
    \item{cohortindex}{\code{signature(x = "APCList")}: 
        get the \code{cohortindex} used in the INLA formula }
    \item{agelabels}{\code{signature(x = "APCList")}: 
        get the \code{agelabels}}
    \item{periodlabels}{\code{signature(x = "APCList")}: 
        get the \code{periodlabels}}
    \item{cohortlabels}{\code{signature(x = "APCList")}: 
        get the \code{cohortlabels}}
    \item{agespec.rate}{\code{signature(x = "APCList")}: 
        get the \code{agespec.rate} slot }
    \item{agespec.proj}{\code{signature(x = "APCList")}: 
        get the \code{agespec.proj} slot }
    \item{agestd.rate}{\code{signature(x = "APCList")}: 
        get the \code{agestd.rate} slot }
    \item{agestd.proj}{\code{signature(x = "APCList")}: 
        get the \code{agestd.proj} slot }
    \item{inlares}{\code{signature(x = "APCList")}: 
        get the \code{inlares} slot }
    \item{periodindex<-}{\code{signature(x = "APCList")}: 
        replace the \code{periodindex} slot used in the INLA formula.
        Caution: Changing this index you should know what you do, 
        otherwise it will result in errors.}
    \item{cohortindex<-}{\code{signature(x = "APCList")}: 
        replace the \code{cohortindex} slot used in the INLA formula.
        Caution: Changing this index you should know what you do, 
        otherwise it will result in errors.}        
    \item{overdisindex<-}{\code{signature(x = "APCList")}: 
        replace the \code{overdisindex} slot used in the INLA formula.
        Caution: Changing this index you should know what you do, 
        otherwise it will result in errors.}               
    \item{agelabels<-}{\code{signature(x = "APCList")}: 
        replace the \code{agelabels} slot }
    \item{periodlabels<-}{\code{signature(x = "APCList")}: 
        replace the \code{periodlabels} slot }
    \item{cohortlabels<-}{\code{signature(x = "APCList")}: 
        replace the \code{cohortlabels} slot }
    \item{agespec.rate<-}{\code{signature(x = "APCList")}: 
        replace the \code{agespec.rate} slot }
    \item{agespec.proj<-}{\code{signature(x = "APCList")}: 
        replace the \code{agespec.proj} slot }
    \item{agestd.rate<-}{\code{signature(x = "APCList")}: 
        replace the \code{agestd.rate} slot }
    \item{agestd.proj<-}{\code{signature(x = "APCList")}: 
        replace the \code{agestd.proj} slot }
    \item{inlares<-}{\code{signature(x = "APCList")}: 
        replace the \code{inlares} slot }
}}
\author{
    Andrea Riebler
}

\examples{
   data(FemLCSweden)
   data(FemPYSweden)

   lc_sweden <- APCList(FemLCSweden, FemPYSweden, gf=5)
   lc_sweden
}
