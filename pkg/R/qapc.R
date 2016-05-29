qapc <- function(APCList, percentiles=c(0.025, seq(0.05, 0.95, by=0.05), 0.975)){
  
  if(class(APCList) != "APCList"){
    stop("The first function argument must be of class \"APCList\"")
  }
  if(!is.vector(percentiles)|| any(!is.numeric(percentiles) | percentiles > 1.0 | percentiles < 0.0)){
    stop("The second function argument must be a numeric vector with values between 0 and 1.")
  }
  
  # we always have age-specific projections
  ap <- agespec.proj(APCList)
  ar <- agespec.rate(APCList)
  
  newcolnames <-  paste(percentiles, "Q", sep="")
  oldcolnames <- colnames(ap[[1]])
  lold <- length(oldcolnames)
  which.new <- which(!(newcolnames %in% oldcolnames))
  if(length(which.new) == 0){
     stop("The desired percentiles are already available.")
  }
  if(any(newcolnames %in% oldcolnames)){
    print("Note: Some or all of the percentiles are already available.")
  }
  percentiles <- percentiles[which.new]
  
  I <- nage(APCList)
  J <- nperiod(APCList)
  res <- inlares(APCList)
  # compute certain quantiles for the linear predictor
  inlaq <- t(sapply(1:(I*J), 
        function(i){INLA::inla.qmarginal(percentiles, res$marginals.linear.predictor[[i]])}))
  # transform to the observation scale
  pre <- exp(inlaq)
  
  for(j in 1:nage(APCList)){
    ap[[j]] = cbind(ap[[j]], 
                    t(sapply(1:J, 
                      function(i){qnorm(percentiles, mean = ap[[j]][i,1], sd=ap[[j]][i,2])})))
    ap[[j]][ap[[j]] < 0] = 0
    if(lold == 2){
      colnames(ap[[j]]) <- c("mean", "sd", paste(percentiles, "Q", sep=""))
    } else {
      colnames(ap[[j]]) <- c("mean", "sd", oldcolnames[-c(1,2)], paste(percentiles, "Q", sep=""))
    }
  }
  agespec.proj(APCList) <- ap
  
  plab=periodlabels(APCList)
  agespec.rate(APCList) <- lapply(1:I, function(m){tmp = cbind(ar[[m]], pre[((m-1)*J+1):(m*J),])
    rownames(tmp) <- plab
    if(lold == 2){
      colnames(tmp) <- c("mean", "sd", paste(percentiles, "Q", sep=""))
    } else {
      colnames(tmp) <- c("mean", "sd", oldcolnames[-c(1,2)], paste(percentiles, "Q", sep=""))
    }
    return(tmp)})
  
  
  ### if we have age-standardized projections
  if(!any(is.na(agestd.proj(APCList)))){
    ## get estimates for mean and sd 
    astdp <- agestd.proj(APCList)
    astdr <- agestd.rate(APCList)
    
    my.quant <- t(sapply(1:J, function(j){qnorm(percentiles, mean=astdr[j,1], sd=astdr[j,2])}))
 
    tmpr = cbind(astdr, my.quant)
    rownames(tmpr) = plab
    
    my.quantp <-  t(sapply(1:J, function(j){qnorm(percentiles, mean=astdp[j,1], sd=astdp[j,2])}))
    tmpp = cbind(astdp, my.quantp)
    rownames(tmpp) = plab
    
    if(lold == 2){
      colnames(tmpr) = c("mean", "sd", paste(percentiles, "Q", sep=""))
      colnames(tmpp) = c("mean", "sd", paste(percentiles, "Q", sep=""))
    } else {
      colnames(tmpr) = c("mean", "sd", oldcolnames[-c(1,2)], paste(percentiles, "Q", sep=""))
      colnames(tmpp) = c("mean", "sd", oldcolnames[-c(1,2)], paste(percentiles, "Q", sep=""))
    }
    agestd.rate(APCList) <- tmpr
    agestd.proj(APCList) <- tmpp
  }
  
  return(APCList)
}

