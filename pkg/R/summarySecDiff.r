summarySecDiff = function(APCList, variable="age", log=FALSE){
  
  result = APCList@inlares[[1]]
 
  if(length(grep("diff_i", rownames(result$summary.lincomb.derived)))==0){
    stop("\n\n\tMarginals for the second differences have not been computed yet.\n\tYou need to switch on the option secondDiff in the BAPC function call.\n")
  }
  if(!(variable %in% c("age", "period", "cohort"))){
    stop("\n\n\tThe varible indicator must be either equal to \"age\", \"period\" or \"cohort\"!\n")
  }
  if(variable == "age"){
    len = nage(APCList)
    idx = grep("diff_i", rownames(result$summary.lincomb.derived))
  } else if(variable == "period"){
    len = nperiod(APCList)
    idx = grep("diff_j", rownames(result$summary.lincomb.derived))
  } else{
    len = ncohort(APCList)
    idx = grep("diff_k", rownames(result$summary.lincomb.derived))
  }
  
  cnames = c("mean", "sd", "0.025Q", "0.5Q", "0.975Q")
  sumTable =  matrix(NA, nrow=len-2, ncol=5)
  colnames(sumTable) = cnames
  
  if(!log){
    for(i in 1:(len-2)){
      sumTable[i,] = unlist(.logsum2expsum(result$marginals.lincomb.derived[[idx[i]]]))[-c(4,6)]
     }
  } else {
    sumTable = result$summary.lincomb.derived[idx,]
  }

  return(sumTable)
}