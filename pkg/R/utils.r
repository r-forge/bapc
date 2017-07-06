.cor2cov <- function(corr_matrix, sd_vector){

    #diag(corr_matrix) <- sd_vector^2
    nc <- ncol(corr_matrix)
    nr <- nrow(corr_matrix)
    cov_matrix <- matrix(NA, nrow=nr, ncol=nc)
    for(i in 1:nr){
        for(j in 1:nc){
            cov_matrix[i,j] <- corr_matrix[i,j] * sd_vector[i] * sd_vector[j]
        }
    }
    return(cov_matrix)
}

.transform_linearPredictor <- function(inla_object){

  # empty list to keep the exp-transformed marginal of the linear predictor
  transformed_marg <- list()
  # get the number of observations
  N <- length(inla_object$.args$data$y)
  # go through all observations
  for(idx in 1:N) {
	# check whether the current observation is NA (to be predicted)
	if(is.na(inla_object$.args$data$y[idx])) {
	  # transform the marginal using the inverse link function, i.e. exp for Poisson
	  # comment: use many points (variable between age groups) to get smooth 
	  # predictive distributions 
	  # (otherwise cyclical patterns will appear  in the tails)
	  ## CAUTION:
	  ## In newer INLA versions method="linear" has to be set explicitly
	  ## to avoid cyclical patterns in the marginals:
	  n=nrow(inla_object$marginals.linear.predictor[[idx]])
	  tmp <- INLA::inla.tmarginal(function(x){exp(x)}, 
	  	inla_object$marginals.linear.predictor[[idx]], method="linear", n=20*n)
	  ## tmp is already a matrix,
	  transformed_marg[[idx]] <- tmp
	} else {
	  # the argument marginals.fitted.values contains the transformed marginals
	  # for non-missing observations
	  transformed_marg[[idx]]<- INLA::inla.smarginal(inla_object$marginals.fitted.values[[idx]], keep.type=TRUE)
	}
  }
  return(transformed_marg)
}


.get_predictions <- function(inla_object, agegroup){

  # get the population counts for the agegroup considered
  pop <- inla_object$.args$data$n[inla_object$.args$data$i==agegroup]
  # just used to get idea of the range of the predictive distribution
  fitted_values <- exp(inla_object$summary.linear.predictor[,c("0.025quant", "0.5quant", "0.975quant")])
  ubound <- floor(fitted_values[inla_object$.args$data$i==agegroup,3]*pop)

  # transform the marginals of the linear predictor to user-scale
  marginals <- .transform_linearPredictor(inla_object)
  # get the indices for the age group considered
  s <- which(inla_object$.args$data$i==agegroup)
  q <- c()
  count <- 1
  for(m in s){
    # just produce some indentifier for the graphics of the whole posterior distribution
    #info <- paste(objectname,": Agegroup =", agegroup, ", index = ", m, "count = ", count)
    # filename to store the whole predictive distribution
    #file <- paste(objectname,"_agegroup_", agegroup, "_index_", m, "_count_", count, ".txt", sep="")
    #print(info)
    # get the right marginal (c is a matrix)
    c <- marginals[[m]]
    # determine the predictive distribution and therefrom the quantiles
    #d <- .predict_dist(pop[count], c, ubound[count], info,file)
    d <- .predict_dist(pop[count], c, ubound[count])
    #print(d)
    # append the quantiles to a matrix
    q <- rbind(q,unlist(d))
    count <- count +1 
  }
  return(q)
}

.predict_dist <- function(pop, param, ubound){
   
  # get the range of plausible y-values
  lu <- 0
  y_range <- lu:max(50, ceiling(2.3*ubound))
  s <- rep(0, length=length(y_range))
 
  # number of evaluation points we have for the marginal
  num_param <- nrow(param)
  num_lambda <- num_param-1

  trapz_vec <- rep(NA, num_lambda)
  lambda_vec <- pop*((param[1:(num_param-1),1] + param[2:num_param,1])/2)
  # trapz_vec necessary to adjust for continuous linear predictor for which
  # we only have finite evaluation points.
  for(k in 1:num_lambda){
	  trapz_vec[k] <- trapz(param[k:(k+1),1], param[k:(k+1),2])
  }

  # for each possible value y_ij
  for(j in 1:length(y_range)){
    # determine the density at point y_range[j], with parameter lambda_vec
    tmp <- dpois(y_range[j], lambda=lambda_vec)*trapz_vec
    # sum the values up 
    s[j] <- sum(tmp) 
  }  

  # The quantile is left continuous: q is the largest integer x such that P(X <= x) < q. 
  quant_sum <- cumsum(s)

  # test whether we got the whole marginal (area should be ~1)
  area <- quant_sum[length(quant_sum)]
  # generate a plot of the posterior density to visually check if everything is ok
#   if(0){
# 	plot(y_range,s, main=paste("Area:", area), sub=info, type="l")
# 	if(area < 0.99 | area > 1.01){
# 	  print(paste("WARNING: Area under the marginal is not 1, but", area))
# 	}
#   }
  # keep the whole predictive distribution
#   if(0){
# 	  write.table(cbind(y_range,s), file=paste("./predDist/", file, sep=""), quote=F, row.names=F, col.names=F)
#   }
  # get the position of the quantiles 0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975 in quant_sum
  quantiles_location <- findInterval(c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975), quant_sum) 
  if(quantiles_location[length(quantiles_location)] == length(quant_sum)){
	  warning("Possibly not enough info: check!!!!!")
  }  
  if(is.element(0, quantiles_location)){
    # at zero the density is larger than 0.025 => location index is zero but there is no
    # element zero and observing 0 counts is the minimum => set the quantile to 0
    quantiles_location[which(quantiles_location==0)] <- 1
  }
  # derive the quantiles
  quantiles <- y_range[quantiles_location]
  names(quantiles) <-  paste(c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975), "quant", sep="")

#   plot(y_range,s, xlim=c(0, max(y_range)), ylim=c(0, 0.2), type="l")
#   cp = inla.tmarginal(function(x){x*pop}, c)
#   lines(cp, col=2)
#   abline(v=quantiles)
  
  return(quantiles)
}

.checkmodel <- function(modlist){

    if(is.null(modlist$include))
        modlist$include = TRUE
    if(is.null(modlist$model))
        modlist$model = "rw2"
    if(is.null(modlist$initial))
        modlist$initial = 4
    if(is.null(modlist$param))
        modlist$param = c(1, 0.00005)
    if(is.null(modlist$prior))
        modlist$prior = "loggamma"
    if(is.null(modlist$scale.model))
        modlist$scale.model = FALSE
    return(modlist)
}

.checkoverdis <- function(modlist){

   if(is.null(modlist$include))
        modlist$include = TRUE
    if(is.null(modlist$model))
        modlist$model = "iid"
    if(is.null(modlist$initial))
        modlist$initial = 4
    if(is.null(modlist$param))
        modlist$param = c(1, 0.00005)
    if(is.null(modlist$prior))
        modlist$prior = "loggamma"
    return(modlist)
}

.num <- function(x, width = if (length(x) > 1) .numlen(x) else 8, 
                          digits = max(4, width)) 
{
  return(formatC(x, format = "g", width = width, flag = "0", 
                 digits = digits))
}

.numlen <- function (n) 
{
  return(floor(log10(max(abs(n)))) + 1)
}

.logsum2expsum = function(marginal){
  
  marginal = INLA::inla.tmarginal(function(x) exp(x), marginal)
  return(INLA::inla.zmarginal(marginal, silent=TRUE))
}

.precsum2varsum = function(marginal){
  
  marginal = INLA::inla.tmarginal(function(x){1/x}, marginal)
  return(INLA::inla.zmarginal(marginal, silent=TRUE))
}
