BAPC <- function(APCList,  predict = list(npredict = 0, retro = TRUE), 
                      model=list(age=list(model="rw2", prior="loggamma", param=c(1, 0.00005), initial=4, scale.model=FALSE), 
                      period=list(include=TRUE, model="rw2", prior="loggamma", param=c(1, 0.00005), initial=4, scale.model=FALSE),
                      cohort=list(include=TRUE, model="rw2", prior="loggamma", param=c(1, 0.00005), initial=4, scale.model=FALSE),
                      overdis=list(include=TRUE, model="iid", prior="loggamma", param=c(1, 0.005), initial=4)),
                      verbose=FALSE, stdweight=NULL){

  if(!is.null(stdweight)){
    if(length(stdweight) != nage(APCList))
        stop("\n\n\tThe weight vector for age-standardization is not of correct length!\n")
  }
  if(!is.null(stdweight) & (sum(stdweight)!= 1)){
    stdweight=stdweight/sum(stdweight)
    cat("\n\n\tCaution: Your age-specific weights were normalised so that they sum to one.\n")
  }  
  if(predict$npredict < 0){
    stop("\n\n\tThe argument npredict must be an integer larger or equal to zero\n")
  }
  model$age <- .checkmodel(model$age)
  model$period <- .checkmodel(model$period)
  model$cohort <- .checkmodel(model$cohort)
  model$overdis <- .checkoverdis(model$overdis)

  APCList@npred <- predict$npredict

  # number of a  data.inla <- data.frame(y=y, n=n, i=i, j=j, k=k, z=z)
  I <- nage(APCList)
  # number of periods
  J <- nperiod(APCList)
  # number of cohorts
  K <- ncohort(APCList)

  y <- APCList@epi
  # retrospective prediction
  if(predict$npredict > 0){
     if(predict$retro){
          y <- epi(APCList)
          y[(J-predict$npredict + 1):J, ] <- NA
          y <- c(y)
    }
  } 
  n <- APCList@pyrs
  # age index
  i <- ageindex(APCList)
  # period index
  j <- periodindex(APCList)
  # cohort index
  k <- cohortindex(APCList)
  # overdispersion parameter
  z <- overdisindex(APCList)

  # data frame for inla
  data.inla <- data.frame(y=y, n=n, i=i, j=j, k=k, z=z)

  # formula object
  formula <- paste("y ~ f(i, model=\"", model$age$model , "\", 
        hyper=list(prec=list(prior=\"", model$age$prior, "\", param=c(", model$age$param[1], ",", model$age$param[2], "), 
            initial=", model$age$initial, ")), constr=T, scale.model=", model$age$scale.model, ")", sep="")

  if(model$period$include){
    if(model$period$model != "drift"){
    formula <- paste(formula, paste("f(j, model=\"", model$period$model , "\", 
        hyper=list(prec=list(prior=\"", model$period$prior, "\", param=c(", model$period$param[1], ",", model$period$param[2], "), 
        initial=", model$period$initial, ")), scale.model=", model$period$scale.model, ")", sep=""), collapse="+", sep="+")
    } else {
      formula <- paste(formula, "j", collapse="+", sep="+")
    }
  }
  if(model$cohort$include){
    if(model$cohort$model != "drift"){
    formula <- paste(formula, paste("f(k, model=\"", model$cohort$model , "\", 
        hyper=list(prec=list(prior=\"", model$cohort$prior, "\", param=c(", model$cohort$param[1], ",", model$cohort$param[2], "), 
        initial=", model$cohort$initial, ")), scale.model=", model$cohort$scale.model, ")", sep=""), collapse="+", sep="+")
    } else {
        formula <- paste(formula, "k", collapse="+", sep="+")
    }
  }
  if(model$overdis$include){
        formula <- paste(formula, paste("f(z, model=\"", model$overdis$model , "\", 
        hyper=list(prec=list(prior=\"", model$overdis$prior, "\", param=c(", model$overdis$param[1], ",", model$overdis$param[2], "), 
        initial=", model$overdis$initial, ")))", sep=""), collapse="+", sep="+")
  }
  if(!is.null(stdweight)){
    # make linear combinations which are the nPred linear predictors
    len <- length(y)
    p <- matrix(NA, len, len)
    diag(p) <- 1
    lcs <- INLA::inla.make.lincombs(Predictor=p)
    names(lcs) <- unlist(strsplit(sprintf("lc%.6f", 1:len/10000), ".", fixed=T))[seq(2, 2*len, 2)]
    config <- TRUE
  } else {
    lcs <- NULL
    config <- TRUE
  }

  res <- INLA::inla(as.formula(formula), family="Poisson", data=data.inla, E=n, 
      control.predictor=list(compute=TRUE),
      lincomb=lcs,
      quantiles=c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975),
      control.compute=list(config=config, dic=TRUE, cpo=TRUE),
      control.inla=list(lincomb.derived.only=TRUE, lincomb.derived.correlation.matrix=TRUE), verbose=verbose)
  inlares(APCList) <- res

  lambda=sapply(1:(I*J), function(i){INLA::inla.emarginal(function(x){exp(x)}, res$marginals.linear.predictor[[i]])})
  lambda2=sapply(1:(I*J), function(i){INLA::inla.emarginal(function(x){exp(x)^2}, res$marginals.linear.predictor[[i]])})

  # variance of lambda
  slambda = lambda2 - lambda^2

  # mean of the predictive distribution
  mu = n*lambda
  # variance of the predictive distribution
  sig = mu + n^2*slambda

  mu_m = matrix(mu, ncol=I, nrow=J, byrow=F)
  lci_m = matrix(mu-1.96*sqrt(sig), ncol=I, nrow=J, byrow=F)
  # do not allow for negative projection rates and set them to zero instead
  lci_m[lci_m < 0] = 0
  uci_m = matrix(mu+1.96*sqrt(sig), ncol=I, nrow=J, byrow=F)
  sd_m = matrix(sqrt(sig), ncol=I, nrow=J, byrow=F)

  pre <- exp(res$summary.linear.predictor[, c("0.025quant", "0.975quant")])
  pre <- cbind("0.025Q"=pre[,1], "mean"=lambda, "0.975Q"=pre[,2], "sd"=sqrt(slambda))

  plab=periodlabels(APCList)
  agespec.rate(APCList) <- lapply(1:I, function(m){tmp=pre[((m-1)*J+1):(m*J),]
      rownames(tmp)=plab
      return(tmp)})
  names(agespec.rate(APCList)) <- agelabels(APCList)
  
  agespec.proj(APCList) <- lapply(1:I, function(m){tmp=cbind("0.025Q"=lci_m[,m], "mean"=mu_m[,m], "0.975Q"=uci_m[,m], "sd"=sd_m[,m]);
     rownames(tmp)=plab
     return(tmp)})
  names(agespec.proj(APCList)) <- agelabels(APCList)
  

  if(!is.null(stdweight)){

    stdweight(APCList) <- stdweight
    my.wm <- matrix(rep(stdweight, each=J), byrow=F, nrow=J)
    stdobs(APCList)=rowSums(my.wm * epi(APCList)/pyrs(APCList), na.rm=T)*rowSums(pyrs(APCList))

    idx.sort <- sort(res$summary.lincomb.derived$ID, index.return=T)$ix
    lc <- res$summary.lincomb.derived[idx.sort,]

    mcor <- res$misc$lincomb.derived.correlation.matrix
    if(sum(colnames(mcor)!=rownames(mcor)) > 0){
        message("WARNING")}

    mcor.idx <- sort(as.numeric(colnames(mcor)), index.return=T)$ix
    mcor <- mcor[mcor.idx,]
    mcor <- mcor[, mcor.idx]

    if(sum(colnames(mcor) != rownames(lc)) > 0){
        message("WARNING")}

    # get the covariance matrix for the linear combinations
    mcov <- .cor2cov(mcor, lc$sd)
    # use the multivariate delta rule to get the covariance matrix
    # for exp(linear combinations)
    D_matrix <- diag(exp(lc$mean))
    mcov_exp <- D_matrix %*% mcov %*% t(D_matrix)

    w_matrix <- matrix(0, ncol=(I*J), nrow=J)
    for(m in 1:J){
        w_matrix[m, j==m] <- stdweight
    }
    ## get the variances of the age standardize quantities
    my_sd <- sqrt(diag(w_matrix %*% mcov_exp %*% t(w_matrix)))

    mi <- res$marginals.lincomb.derived[idx.sort]

    #new_mean <- matrix(sapply(1:(I*J), function(u){inla.emarginal(function(y){exp(y)}, mi[[u]])}), nrow=J, byrow=F)
    new_mean <- exp(lc$mean)
    my.wm <- matrix(rep(stdweight, each=J), byrow=F, nrow=J)
    new_mean <- rowSums(my.wm*new_mean)
   
    new_l <- (new_mean - 1.96*my_sd) 
    new_u <- (new_mean + 1.96*my_sd) 

    tmp = cbind("0.025Q"=new_l, "mean"=new_mean, "0.975Q"=new_u, "sd"=my_sd)
    rownames(tmp) = plab
    agestd.rate(APCList) <- tmp
    
    agg.n <- rowSums(pyrs(APCList))
    agg.mean <- agg.n * new_mean
    agg.std <- sqrt(agg.mean + agg.n^2*my_sd^2)
    tmp = cbind("0.025Q"=agg.mean-1.96*agg.std, "mean"=agg.mean,
        "0.975Q"=agg.mean+1.96*agg.std, "sd"=agg.std)
    rownames(tmp) = plab
    agestd.proj(APCList) <- tmp
  }

  return(APCList)
}

