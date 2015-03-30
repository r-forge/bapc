ratesByAge <- function(APCList, age, per, scale=10^3, type="pa", log.ax="", ann=TRUE, grid=TRUE, ...){

    if(!(type %in% c("pa", "ca")))
        stop("Argument type must be either equal to pa or ca")

    np <- nperiod(APCList)
    rates <- epi(APCList)/pyrs(APCList)

    rates <- t(rates)
    rateplot(rates*scale, age=age, per=per, which=type, log.ax=log.ax, ann=ann, 
        grid=grid, ylab=paste("Rates per", format(scale, scientific=FALSE), sep=" "), xlab="Date of event",...)
}
