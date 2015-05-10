plotBAPC <- function(APCList, scale=10^3, type="ageSpecProj", showdata=FALSE, mfrow=NULL, coladd=NULL, ...){

    if(!(type %in% c("ageSpecRate", "ageSpecProj", "ageSpecBoth", "ageStdRate", "ageStdProj", "ageStdBoth")))
      stop("type must be one of ageSpecRate, ageSpecProj, ageSpecBoth, ageStdRate, ageStdProj or ageStdBoth")

    if(type=="ageSpecBoth"){
      projrate <- agespec.rate(APCList)
      projcounts <- agespec.proj(APCList)
      plab <- periodlabels(APCList)

      I <- nage(APCList)
      J <- nperiod(APCList)

      y <- epi(APCList)
      n <- pyrs(APCList)

      if(is.null(mfrow))
        mfrow <- c(ceiling(I/4),4)

      if(is.null(coladd))
        coladd=2

      par(mfrow=mfrow, ...)
      for(u in 1:I){
        if(showdata){
            y_max <- max(y[,u]/n[,u]*scale, as.matrix(projrate[[u]])*scale, projcounts[[u]]/n[,u]*scale, na.rm=TRUE)
        } else {
            y_max <- max(projrate[[u]]*scale, projcounts[[u]]/n[,u]*scale)
        } 
        al <- agelabels(APCList)
        matplot(plab, projrate[[u]]*scale,
	  main=paste("Age group", al[u]), type="n", lty=c(2,1,2), col=1,
	  xlab="Period", ylab=paste("Rate per ", format(scale, scientific=FALSE), sep=""), ylim=c(0, y_max), ...)
        # polygon(c(plab, rev(plab)), c(projcounts[[u]][,1]/n[,u]*scale, rev(projcounts[[u]][,3]/n[,u]*scale)), col=coladd, border=NA)
        polygon(c(plab, rev(plab)), c(projrate[[u]][,1]*scale, rev(projrate[[u]][,3]*scale)), col=coladd, border=NA)
        matplot(plab, projcounts[[u]][,1:3]/n[,u]*scale, lty=c(2,1,2), col=1, add=T, type="l", ...)
        lines(plab,  projrate[[u]][,2]*scale)
        if(showdata){
            if(is.na(sum(y))){
                matplot(plab[1:(J-APCList@npred)], (y[,u]/n[,u]*scale)[1:(J-APCList@npred)], type="p", pch=20, add=T)
            } else {
                matplot(plab, y[,u]/n[,u]*scale, type="p", pch=20, add=T)
            }
        }
        #abline(v=max(as.numeric(plab))-APCList@npred + .5, col="grey55", lwd=1.2)
        numplab <- as.numeric(plab) 
	abline(v=max(numplab)-(APCList@npred*(numplab[2]-numplab[1])) + .5, col="grey55", lwd=1.2)
      }
    }

    if(type=="ageSpecRate"){
      projrate <- agespec.rate(APCList)
      plab <- periodlabels(APCList)

      I <- nage(APCList)
      J <- nperiod(APCList)

      y <- epi(APCList)
      n <- pyrs(APCList)

      if(is.null(mfrow))
        mfrow <- c(ceiling(I/4),4)

      par(mfrow=mfrow, ...)
      for(u in 1:I){
        if(showdata){
            y_max <- max(y[,u]/n[,u]*scale, as.matrix(projrate[[u]][,1:3])*scale, na.rm=TRUE)
        } else {
            y_max <- max(projrate[[u]][,1:3]*scale)
        } 
        al <- agelabels(APCList)
        matplot(plab, projrate[[u]][,1:3]*scale,
        main=paste("Age group", al[u]), type="l", lty=c(2,1,2), col=1,
        xlab="Period", ylab=paste("Rate per ",format(scale, scientific=FALSE), sep=""), ylim=c(0, y_max), ...)
        if(showdata){
            if(is.na(sum(y))){
                matplot(plab[1:(J-APCList@npred)], (y[,u]/n[,u]*scale)[1:(J-APCList@npred)], type="p", pch=20, add=T)
            } else {
                matplot(plab, y[,u]/n[,u]*scale, type="p", pch=20, add=T)
            }
        }
        #abline(v=max(as.numeric(plab))-APCList@npred + .5, col="grey55", lwd=1.2)
        numplab <- as.numeric(plab) 
	abline(v=max(numplab)-(APCList@npred*(numplab[2]-numplab[1])) + .5, col="grey55", lwd=1.2)

      }
    }
    if(type=="ageSpecProj"){
	  projcounts <- agespec.proj(APCList)
	  plab <- periodlabels(APCList)

	  I <- nage(APCList)
	  J <- nperiod(APCList)

	  y <- epi(APCList)
	  n <- pyrs(APCList)

	  if(is.null(mfrow))
	    mfrow <- c(ceiling(I/4),4)

	  par(mfrow=mfrow, ...)
	  for(u in 1:I){
          if(showdata){
            y_max <- max(y[,u]/n[,u]*scale,projcounts[[u]]/n[,u]*scale, na.rm=TRUE)
          } else {
            y_max <- max(projcounts[[u]]/n[,u]*scale)
          } 
	      al <- agelabels(APCList)
	      matplot(plab, projcounts[[u]][,1:3]/n[,u]*scale,
	      main=paste("Age group", al[u]), type="l", lty=c(2,1,2), col=1,
	      xlab="Period", ylab=paste("Rate per ", format(scale, scientific=FALSE), sep=""), ylim=c(0, y_max), ...)
	      if(showdata){
            if(is.na(sum(y))){
                matplot(plab[1:(J-APCList@npred)], (y[,u]/n[,u]*scale)[1:(J-APCList@npred)], type="p", pch=20, add=T)
            } else {
                matplot(plab, y[,u]/n[,u]*scale, type="p", pch=20, add=T)
            }
	      }
	      
	     #abline(v=max(as.numeric(plab))-APCList@npred + .5, col="grey55", lwd=1.2)
	             numplab <- as.numeric(plab) 
	abline(v=max(numplab)-(APCList@npred*(numplab[2]-numplab[1])) + .5, col="grey55", lwd=1.2)

	  }
    } 
    if(type=="ageStdBoth"){
      J <- nperiod(APCList)
      plab <- periodlabels(APCList)
      my.wm <- matrix(rep(stdweight(APCList), each=J), byrow=F, nrow=J)
      agg.n <- rowSums(pyrs(APCList))

      if(is.null(coladd))
        coladd=2

      matplot(plab, agestd.rate(APCList)*scale, type="l", lty=c(2,1,2), col=1,  
          xlab="Period", ylab=paste("Agestd. rate per ", format(scale, scientific=FALSE), sep=""),  ...)
      polygon(c(plab, rev(plab)), c(agestd.rate(APCList)[,1]*scale, rev(agestd.rate(APCList)[,3]*scale)), col=coladd, border=NA)
      matplot(plab, agestd.proj(APCList)[,1:3]/agg.n*scale, type="l", lty=c(2,1,2), col=c(1,1,1),  
          add=T, ...)
     lines(plab,  agestd.rate(APCList)[,2]*scale)
      if(showdata){
          y <- epi(APCList)
          n <- pyrs(APCList)
          if(is.na(sum(y))){
            matplot(plab[1:(J-APCList@npred)], rowSums(my.wm * y/n*scale)[1:(J-APCList@npred)], type="p", pch=20, add=T)
          } else {
            matplot(plab, rowSums(my.wm * y/n*scale), type="p", pch=20, add=T)
          }
          
          #abline(v=max(as.numeric(plab))-APCList@npred + .5, col="grey55", lwd=1.2)
                  numplab <- as.numeric(plab) 
	abline(v=max(numplab)-(APCList@npred*(numplab[2]-numplab[1])) + .5, col="grey55", lwd=1.2)

      }
    }

    if(type=="ageStdRate"){
	  J <- nperiod(APCList)
	  plab <- periodlabels(APCList)
	  my.wm <- matrix(rep(stdweight(APCList), each=J), byrow=F, nrow=J)

	  matplot(plab, agestd.rate(APCList)[,1:3]*scale, type="l", lty=c(2,1,2), col=1,  
	      xlab="Period", ylab=paste("Agestd. rate per ", format(scale, scientific=FALSE), sep=""), ...)
	  if(showdata){
	      y <- epi(APCList)
	      n <- pyrs(APCList)
          if(is.na(sum(y))){
            matplot(plab[1:(J-APCList@npred)], rowSums(my.wm * y/n*scale)[1:(J-APCList@npred)], type="p", pch=20, add=T)
          } else {
            matplot(plab, rowSums(my.wm * y/n*scale), type="p", pch=20, add=T)
          }
	      #abline(v=max(as.numeric(plab))-APCList@npred + .5, col="grey55", lwd=1.2)
	              numplab <- as.numeric(plab) 
	abline(v=max(numplab)-(APCList@npred*(numplab[2]-numplab[1])) + .5, col="grey55", lwd=1.2)

      }
    }
    if(type=="ageStdProj"){
      J <- nperiod(APCList)
      plab <- periodlabels(APCList)
      my.wm <- matrix(rep(stdweight(APCList), each=J), byrow=F, nrow=J)
      agg.n <- rowSums(pyrs(APCList))
      matplot(plab, agestd.proj(APCList)[,1:3]/agg.n*scale, type="l", lty=c(2,1,2), col=1,  
          xlab="Period", ylab=paste("Agestd. rate per ", format(scale, scientific=FALSE), sep=""), ...)
      if(showdata){
          y <- epi(APCList)
          n <- pyrs(APCList)
          if(is.na(sum(y))){
            matplot(plab[1:(J-APCList@npred)], rowSums(my.wm * y/n*scale, na.rm=T)[1:(J-APCList@npred)], type="p", pch=20, add=T)
          } else {
            matplot(plab, rowSums(my.wm * y/n*scale, na.rm=T), type="p", pch=20, add=T)
          }
         #abline(v=max(as.numeric(plab))-APCList@npred + .5, col="grey55", lwd=1.2)
                  numplab <- as.numeric(plab) 
	abline(v=max(numplab)-(APCList@npred*(numplab[2]-numplab[1])) + .5, col="grey55", lwd=1.2)

      }
    }
}

