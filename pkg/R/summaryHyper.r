summaryHyper = function(APCList, var=TRUE){
  
  result = APCList@inlares[[1]]
  dim = nrow(result$summary.hyperpar)
  if(var){
    sumTable = matrix(NA, nrow=dim, ncol=5)
  colnames(sumTable) = c("mean", "sd", "0.025Q", "0.5Q", "0.975Q")
    precnames = rownames(result$summary.hyperpar)
    rownames(sumTable) = gsub("Precision","Variance", precnames)
    for(i in 1:dim){
      sumTable[i,] = unlist(.precsum2varsum(result$marginals.hyperpar[[i]]))[-c(4,6)]
    }
  } else {
      sumTable = result$summary.hyperpar
  }
  return(sumTable)
}