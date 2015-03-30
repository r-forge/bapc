setClass("APCList", representation(
    epi="numeric", 
    pyrs="numeric",  
    gf="numeric",
    agelab="character",
    periodlab="character",
    cohortlab="character",
    age.idx="numeric",
    period.idx="numeric",
    cohort.idx="numeric",
    stdweight="numeric",
    stdobs="numeric",
    npred="numeric",
    agespec.rate="list",
    agespec.proj="list",
    agestd.rate="matrix",
    agestd.proj="matrix",
    inlares="list"
    ))

setGeneric("APCList", function(epi, pyrs, gf, ... )
           {standardGeneric("APCList")})

setMethod("APCList", c("data.frame", "data.frame", "numeric"),
    function(epi, pyrs, gf, agelab=character(), periodlab=character(), cohortlab=character(),
    stdweight=numeric(), stdobs=numeric(), npred=numeric(), agespec.rate=list(), agespec.proj=list(), 
    agestd.rate=matrix(), agestd.proj=matrix())
{
    if((class(epi)!= "data.frame") || (class(pyrs) != "data.frame") || (class(gf) != "numeric")){
            stop("\n\n\t  `epi' and `pyrs' must be of class `data.frame' and `gf' a `numeric' number\n\n")
    }
    if(is.null(epi) || is.null(pyrs) || is.null(gf)){
            stop("\n\n\t  `epi' and `pyrs' and `gf' are obligatory and must be provided \n\n")
    }
    
    if(any(is.na(pyrs))){
       stop("\n\n\t There are no NAs in pyrs allowed\n\n")
    }
    ## get the object dimensions
    nad <- ncol(epi)
    npd <- nrow(epi)
    nap <- ncol(pyrs)
    npp <- nrow(pyrs)

    if((nad != nap) || (npd != npp))
        stop("\n\n\tepi and pyrs must be of same dimension!\n\n")
    na <- length(agelab)
    np <- length(periodlab)
    nc <- length(cohortlab)

    if((na > 0) && (na != nad))
        stop("\n\n\tage labels have not the correct dimension!\n\n")
    if((np > 0) && (np != npd))
        stop("\n\n\tperiod labels have not the correct dimension!\n\n")
    if((nc > 0) && (na != (gf*(na-1)+np)))
        stop("\n\n\tcohort labels have not the correct dimension!\n\n")
    if(na == 0)
        agelab <- colnames(epi)
    if(length(stdweight) > 0 && length(stdweight) != nad)
        stop("\n\n\tthe vector with the standardized weights has not the correct dimension\n\n")
   
    if(np == 0)
        periodlab <- rownames(epi)

    if(nc == 0)
        cohortlab <- as.character(1:(gf*(nap - 1)+npp))


    epi <- c(as.matrix(epi))
    pyrs <- c(as.matrix(pyrs))
    age.idx <- rep(1:nap, each=npp)
    period.idx <- rep(1:npp, nap)
    cohort.idx <- gf*(nap - age.idx) + period.idx
    stdweight <- stdweight/sum(stdweight)

    new("APCList", epi=epi, pyrs=pyrs, gf=gf, 
        agelab=agelab, periodlab=periodlab, cohortlab=cohortlab,
        age.idx=age.idx, period.idx=period.idx, cohort.idx=cohort.idx,
        stdweight=stdweight, stdobs=stdobs, npred=npred, agespec.rate=agespec.rate, 
        agespec.proj=agespec.proj, agestd.rate=agestd.rate, 
        agestd.proj=agestd.proj)
})

#######################################################
## Show function for the "APCList" class
######################################################################

setMethod("show", "APCList", function(object) {
    cat("Object of class 'APCList'.\n\n")
    cat("- Number of age groups:", length(object@agelab), "\n\n")
    cat("- Number of periods:", length(object@periodlab), "\n\n")
    cat("- Gridfactor:", object@gf, "\n\n")
})

######################################################################
## Access functions for the "APCList" class
######################################################################

if(!isGeneric("[")) setGeneric("[", function(object) standardGeneric("["))
setMethod("[", signature(x="APCList",i="ANY",j="ANY"),
    function(x, i, j, ..., drop=TRUE) {
    
    y <- matrix(x@epi, nrow=length(x@periodlab), byrow=F)
    n <- matrix(x@pyrs, nrow=length(x@periodlab), byrow=F)

    y <- y[i, j]
    n <- n[i, j]

    nap <- ncol(y)
    npp <- nrow(y)

    agelab <- x@agelab[j]
    periodlab <- x@periodlab[i]
    cohortlab <- as.character(1:(x@gf*(nap- 1)+npp))

    age.idx <- rep(1:nap, each=npp)
    period.idx <- rep(1:npp, nap)
    cohort.idx <- x@gf*(nap - age.idx) + period.idx

    message("\n\tAll projections were deleted\n\n")
    new("APCList", epi=c(y), pyrs=c(n), gf=x@gf, 
        agelab=agelab, periodlab=periodlab, cohortlab=cohortlab,
        age.idx=age.idx, period.idx=period.idx, cohort.idx=cohort.idx,
        stdweight=numeric(), stdobs=numeric(), npred=numeric(), 
        agespec.rate=list(), agespec.proj=list(), 
        agestd.rate=matrix(), agestd.proj=matrix())
})

if(!isGeneric("nage")) setGeneric("nage", 
    function(x) standardGeneric("nage"))
setMethod("nage", "APCList", function(x) {
        length(x@agelab)
})
if(!isGeneric("nperiod")) setGeneric("nperiod", 
    function(x) standardGeneric("nperiod"))
setMethod("nperiod", "APCList",
    function(x) {
        length(x@periodlab)
})
if(!isGeneric("ncohort")) setGeneric("ncohort", 
    function(x) standardGeneric("ncohort"))
setMethod("ncohort", "APCList",
    function(x) {
        length(x@cohortlab)
})
if(!isGeneric("gridfactor")) setGeneric("gridfactor", 
    function(x) standardGeneric("gridfactor"))
setMethod("gridfactor", "APCList",
    function(x) {
        x@gf
})
if(!isGeneric("epi")) setGeneric("epi", 
    function(x) standardGeneric("epi"))
setMethod("epi", "APCList",
    function(x) {
        tmp=x@epi
        tmp=matrix(tmp,ncol=nage(x), byrow=FALSE)
        colnames(tmp)=x@agelab
        rownames(tmp)=x@periodlab
        return(tmp)
})
if(!isGeneric("pyrs")) setGeneric("pyrs", 
    function(x) standardGeneric("pyrs"))
setMethod("pyrs", "APCList",
    function(x) {
        tmp=x@pyrs
        tmp=matrix(tmp,ncol=nage(x), byrow=FALSE)
        colnames(tmp)=x@agelab
        rownames(tmp)=x@periodlab
        return(tmp)
})
if(!isGeneric("agelabels")) setGeneric("agelabels", 
    function(x) standardGeneric("agelabels"))
setMethod("agelabels", "APCList",
    function(x) {
        x@agelab
})
if(!isGeneric("periodlabels")) setGeneric("periodlabels", 
    function(x) standardGeneric("periodlabels"))
setMethod("periodlabels", "APCList",
    function(x) {
        x@periodlab
})
if(!isGeneric("cohortlabels")) setGeneric("cohortlabels", 
    function(x) standardGeneric("cohortlabels"))
setMethod("cohortlabels", "APCList",
    function(x) {
        x@cohortlab
})
if(!isGeneric("ageindex")) setGeneric("ageindex", 
    function(x) standardGeneric("ageindex"))
setMethod("ageindex", "APCList",
    function(x) {
        x@age.idx
})
if(!isGeneric("periodindex")) setGeneric("periodindex", 
    function(x) standardGeneric("periodindex"))
setMethod("periodindex", "APCList",
    function(x) {
        x@period.idx
})
if(!isGeneric("cohortindex")) setGeneric("cohortindex", 
    function(x) standardGeneric("cohortindex"))
setMethod("cohortindex", "APCList",
    function(x) {
        x@cohort.idx
})
if(!isGeneric("stdweight")) setGeneric("stdweight", 
    function(x) standardGeneric("stdweight"))
setMethod("stdweight", "APCList",
    function(x) {
        x@stdweight
})
if(!isGeneric("stdobs")) setGeneric("stdobs", 
    function(x) standardGeneric("stdobs"))
setMethod("stdobs", "APCList",
    function(x) {
        x@stdobs
})
if(!isGeneric("agespec.rate")) setGeneric("agespec.rate", 
    function(x) standardGeneric("agespec.rate"))
setMethod("agespec.rate", "APCList",
    function(x) {
        x@agespec.rate
})
if(!isGeneric("agespec.proj")) setGeneric("agespec.proj", 
    function(x) standardGeneric("agespec.proj"))
setMethod("agespec.proj", "APCList",
    function(x) {
        x@agespec.proj
})
if(!isGeneric("agestd.rate")) setGeneric("agestd.rate", 
    function(x) standardGeneric("agestd.rate"))
setMethod("agestd.rate", "APCList",
    function(x) {
        x@agestd.rate
})
if(!isGeneric("agestd.proj")) setGeneric("agestd.proj", 
    function(x) standardGeneric("agestd.proj"))
setMethod("agestd.proj", "APCList",
    function(x) {
        x@agestd.proj
})
if(!isGeneric("inlares")) setGeneric("inlares", 
    function(x) standardGeneric("inlares"))
setMethod("inlares", "APCList",
    function(x) {
        x@inlares[[1]]
})

setGeneric("periodindex<-", function(x, value) standardGeneric("periodindex<-"))
setReplaceMethod("periodindex", "APCList", function(x, value) {
    if(length(value) != length(x@period.idx))
        stop("\n\n\tThe index vector does not have the correct length!\n")
    x@period.idx <- value
    x
})
setGeneric("cohortindex<-", function(x, value) standardGeneric("cohortindex<-"))
setReplaceMethod("cohortindex", "APCList", function(x, value) {
    if(length(value) != length(x@cohort.idx))
        stop("\n\n\tThe index vector does not have the correct length!\n")
    x@cohort.idx <- value
    x
})
setGeneric("agelabels<-", function(x, value) standardGeneric("agelabels<-"))
setReplaceMethod("agelabels", "APCList", function(x, value) {
    if(length(value) != length(x@agelab))
        stop("\n\n\tThe label vector does not have the correct length!\n")
    x@agelab <- value
    x
})
setGeneric("periodlabels<-", function(x, value) standardGeneric("periodlabels<-"))
setReplaceMethod("periodlabels", "APCList", function(x, value) {
    if(length(value) != length(x@periodlab))
        stop("\n\n\tThe label vector does not have the correct length!\n")

    x@periodlab <- value
    x
})
setGeneric("cohortlabels<-", function(x, value) standardGeneric("cohortlabels<-"))
setReplaceMethod("cohortlabels", "APCList", function(x, value) {
    if(length(value) != length(x@cohortlab))
        stop("\n\n\tThe label vector does not have the correct length!\n")
    x@cohortlab <- value
    x
})
setGeneric("stdweight<-", function(x, value) standardGeneric("stdweight<-"))
setReplaceMethod("stdweight", "APCList", function(x, value) {
    x@stdweight <- value
    x
})
setGeneric("stdobs<-", function(x, value) standardGeneric("stdobs<-"))
setReplaceMethod("stdobs", "APCList", function(x, value) {
    x@stdobs <- value
    x
})
setGeneric("agespec.rate<-", function(x, value) standardGeneric("agespec.rate<-"))
setReplaceMethod("agespec.rate", "APCList", function(x, value) {
    x@agespec.rate <- value
    x
})

setGeneric("agespec.proj<-", function(x, value) standardGeneric("agespec.proj<-"))
setReplaceMethod("agespec.proj", "APCList", function(x, value) {
    x@agespec.proj <- value
    x
})

setGeneric("agestd.rate<-", function(x, value) standardGeneric("agestd.rate<-"))
setReplaceMethod("agestd.rate", "APCList", function(x, value) {
    x@agestd.rate <- value
    x
})


setGeneric("agestd.proj<-", function(x, value) standardGeneric("agestd.proj<-"))
setReplaceMethod("agestd.proj", "APCList", function(x, value) {
    x@agestd.proj <- value
    x
})

setGeneric("inlares<-", function(x, value) standardGeneric("inlares<-"))
setReplaceMethod("inlares", "APCList", function(x, value) {
    x@inlares <- list(value)
    x
})


