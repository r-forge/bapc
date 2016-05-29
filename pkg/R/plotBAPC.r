plotBAPC <- function(APCList, scale=10^5, type="ageSpecProj", probs=c(0.025, seq(0.05, 0.95, by=0.05), 0.975), 
                     start=NULL, showdata=TRUE, mfrow=NULL, 
                     col.fan=sequential_hcl, obs.pch=21, obs.col="white", obs.bg="black",  obs.lwd=2, obs.cex=1,
                     ln=NULL,...){

  if(class(APCList) != "APCList"){
    stop("The first function argument must be of class \"APCList\"")
  }
  if(!is.vector(probs)|| any(!is.numeric(probs) | probs > 1.0 | probs < 0.0)){
    stop("The second function argument must be a numeric vector with values between 0 and 1.")
  }
  if(!(type %in% c("ageSpecRate", "ageSpecProj", "ageStdRate", "ageStdProj")))
      stop("type must be one of ageSpecRate, ageSpecProj, ageStdRate, ageStdProj")
  
  newcolnames <-  paste(probs, "Q", sep="")
  if(type=="ageSpecRate"){
      APCList <- qapc(APCList, percentiles=probs)
    
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
        if(is.null(start)){
          st <- 1
        } else {
          if(!is.element(start, periodlabels(APCList))){
            stop("start must be an element of the the period labels!")
          }
          st <- which(periodlabels(APCList)==start)
        }
        if(showdata){
            y_max <- max(y[st:J,u]/n[st:J,u]*scale, 
                         as.matrix(projrate[[u]][st:J,newcolnames])*scale, 
                         na.rm=TRUE)
        } else {
            y_max <- max(projrate[[u]][st:J,newcolnames]*scale)
        } 
        al <- agelabels(APCList)
        # do not plot the x-axis
        matplot(projrate[[u]][st:J,newcolnames]*scale,
          main=paste("Age group", al[u]), type="n", lty=c(2,1,2), col=1,
          xlab="Period", ylab=paste("Rate per ",format(scale, scientific=FALSE), sep=""), 
          ylim=c(0, y_max), xaxt="n", ...)
        # find out where to put the x-axis tick marks and label with plab
        plabn <- as.numeric(plab[st:J])
        getlabels <- cbind(1:length(plabn), plabn)
        plabp <- getlabels[getlabels[,2] %in% pretty(range(plabn)),]
        axis(1, at=plabp[,1], labels=plabp[,2])
        # add fanplot
        fan(t(projrate[[u]][st:J,newcolnames])*scale, data.type='values', probs=probs, 
            fan.col=col.fan, ln=ln)
        lines(projrate[[u]][st:J,1]*scale, col=1)
        if(showdata){
            if(is.na(sum(y))){
                matplot((y[,u]/n[,u]*scale)[st:(J-APCList@npred)], type="p", 
                    pch=obs.pch, col=obs.col, bg=obs.bg, lwd=obs.lwd, cex=obs.cex, add=T)
            } else {
                matplot(y[st:J,u]/n[st:J,u]*scale, type="p", 
                    pch=obs.pch, col=obs.col, bg=obs.bg, lwd=obs.lwd, cex=obs.cex, add=T)
            }
        }
        abline(v=length(plab)-APCList@npred - st + 1.5, col="grey55", lty=2, lwd=1.2)
        #numplab <- as.numeric(plab) 
	      #abline(v=max(numplab)-(APCList@npred*(numplab[2]-numplab[1])) + start -.5, col=1, lwd=1.2)
      }
    }
    if(type=="ageSpecProj"){
      APCList <- qapc(APCList, percentiles=probs)
      
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
	        if(is.null(start)){
	         st <- 1
	        } else {
	          if(!is.element(start, periodlabels(APCList))){
	            stop("start must be an element of the the period labels!")
	          }
	          st <- which(periodlabels(APCList)==start)
	        }
	      
          if(showdata){
            if(is.na(sum(y))){
              y_max <- max(y[st:J,u]/n[st:J,u]*scale, projcounts[[u]][st:J,newcolnames]/n[st:J,u]*scale, na.rm=TRUE)
            } else {
              y_max <- max(y[st:J,u]/n[st:J,u]*scale, projcounts[[u]][st:J,newcolnames]/n[st:J,u]*scale)
            }
          } 
	      al <- agelabels(APCList)
	      matplot(projcounts[[u]][st:J,newcolnames]/n[st:J,u]*scale,
	        main=paste("Age group", al[u]), type="n", lty=c(2,1,2), col=1,
	        xlab="Period", ylab=paste("Rate per ", format(scale, scientific=FALSE), sep=""), ylim=c(0, y_max),
	        xaxt="n", ...)
	      # find out where to put the x-axis tick marks and label with plab
	      plabn <- as.numeric(plab[st:J])
	      getlabels <- cbind(1:length(plabn), plabn)
	      plabp <- getlabels[getlabels[,2] %in% pretty(range(plabn)),]
	      axis(1, at=plabp[,1], labels=plabp[,2])
	      # add fanplot
	      fan(t(projcounts[[u]][st:J,newcolnames]/n[st:J,u])*scale, data.type='values', probs=probs, 
	          fan.col=col.fan, ln=ln)
	      lines(projcounts[[u]][st:J,1]/n[st:J,u]*scale, col=1)
	      
	      if(showdata){   
            if(is.na(sum(y))){
                matplot((y[,u]/n[,u]*scale)[st:(J-APCList@npred)], type="p", 
                        pch=obs.pch, col=obs.col, bg=obs.bg, lwd=obs.lwd,cex=obs.cex,  add=T)
            } else {
                matplot(y[st:J,u]/n[st:J,u]*scale,  type="p", 
                        pch=obs.pch, col=obs.col, bg=obs.bg, lwd=obs.lwd,cex=obs.cex,  add=T)
            }
	      }
	      
	     #abline(v=max(as.numeric(plab))-APCList@npred + .5, col="grey55", lwd=1.2)
	     #        numplab <- as.numeric(plab) 
	     abline(v=length(plab)-APCList@npred - st + 1.5, col="grey55", lty=2, lwd=1.2)
	    }
    }
    if(type == "ageStdRate" | type == "ageStdProj"){
      if(is.na(agestd.proj(APCList))){
        stop("There are no agestandardized projections or rates available. Probably you did not provide standardisation weights to the function BAPC\n")
      }
    if(type=="ageStdRate"){
      APCList <- qapc(APCList, percentiles=probs)
      
	    J <- nperiod(APCList)
	    plab <- periodlabels(APCList)
	    my.wm <- matrix(rep(stdweight(APCList), each=J), byrow=F, nrow=J)

	    if(is.null(start)){
	      st <- 1
	    } else {
	      if(!is.element(start, periodlabels(APCList))){
	        stop("start must be an element of the the period labels!")
	      }
	      st <- which(periodlabels(APCList)==start)
	    }
	    matplot(agestd.rate(APCList)[st:J,newcolnames]*scale, type="n", lty=c(2,1,2), col=1,  
	      xlab="Period", ylab=paste("Agestd. rate per ", format(scale, scientific=FALSE), sep=""), xaxt="n",...)
	    # find out where to put the x-axis tick marks and label with plab
	    plabn <- as.numeric(plab[st:J])
	    getlabels <- cbind(1:length(plabn), plabn)
	    plabp <- getlabels[getlabels[,2] %in% pretty(range(plabn)),]
	    axis(1, at=plabp[,1], labels=plabp[,2])
	    # add fanplot
	    fan(t(agestd.rate(APCList)[st:J,newcolnames]*scale), data.type='values', probs=probs, 
	        fan.col=col.fan, ln=ln)
	    lines(agestd.rate(APCList)[st:J,1]*scale, col=1)
	    
	    if(showdata){
	       y <- epi(APCList)
	       n <- pyrs(APCList)
          if(is.na(sum(y))){
            matplot(rowSums(my.wm * y/n*scale)[st:(J-APCList@npred)], 
                    type="p",pch=obs.pch, col=obs.col, bg=obs.bg, lwd=obs.lwd,cex=obs.cex, add=T)
          } else {
            matplot(rowSums(my.wm * y/n*scale)[st:J], type="p", 
                    pch=obs.pch, col=obs.col, bg=obs.bg, lwd=obs.lwd, cex=obs.cex, add=T)
          }
	        #abline(v=max(as.numeric(plab))-APCList@npred + .5, col="grey55", lwd=1.2)
	       #       numplab <- as.numeric(plab) 
	       abline(v=length(plab)-APCList@npred - st + 1.5, col="grey55", lty=2, lwd=1.2)
	    }
    }
    if(type=="ageStdProj"){
      APCList <- qapc(APCList, percentiles=probs)
      
      J <- nperiod(APCList)
      plab <- periodlabels(APCList)
      my.wm <- matrix(rep(stdweight(APCList), each=J), byrow=F, nrow=J)
      agg.n <- rowSums(pyrs(APCList))
      
      if(is.null(start)){
        st <- 1
      } else {
        if(!is.element(start, periodlabels(APCList))){
          stop("start must be an element of the the period labels!")
        }
        st <- which(periodlabels(APCList)==start)
      }
      matplot(agestd.proj(APCList)[st:J,newcolnames]/agg.n[st:J]*scale, type="n", lty=c(2,1,2), 
          col=1,  
          xlab="Period", ylab=paste("Agestd. rate per ", format(scale, scientific=FALSE), sep=""),
          xaxt="n", ...)
      # find out where to put the x-axis tick marks and label with plab)
      plabn <- as.numeric(plab[st:J])
      getlabels <- cbind(1:length(plabn), plabn)
      plabp <- getlabels[getlabels[,2] %in% pretty(range(plabn)),]
      axis(1, at=plabp[,1], labels=plabp[,2])
      # add fanplot
      fan(t(agestd.proj(APCList)[st:J,newcolnames]/agg.n[st:J]*scale), data.type='values', probs=probs, 
          fan.col=col.fan, ln=ln)
      lines(agestd.proj(APCList)[st:J,1]/agg.n[st:J]*scale, col=1)
      
      if(showdata){
          y <- epi(APCList)
          n <- pyrs(APCList)
          if(is.na(sum(y))){
            matplot(rowSums(my.wm * y/n*scale, na.rm=T)[st:(J-APCList@npred)], 
                    pch=obs.pch, col=obs.col, bg=obs.bg, lwd=obs.lwd, cex=obs.cex, type="p", add=T)
          } else {
            matplot(rowSums(my.wm * y/n*scale, na.rm=T)[st:J], 
                    pch=obs.pch, col=obs.col, bg=obs.bg, lwd=obs.lwd, cex=obs.cex, type="p", add=T)
          }
         #abline(v=max(as.numeric(plab))-APCList@npred + .5, col="grey55", lwd=1.2)
         #         numplab <- as.numeric(plab) 
         abline(v=length(plab)-APCList@npred - st + 1.5, col="grey55", lty=2, lwd=1.2)
      }
    }
    }
}
# 
# ## old
# plotBAPC <- function(APCList, scale=10^3, type="ageSpecProj", showdata=FALSE, mfrow=NULL, coladd=NULL, ...){
#   
#   if(!(type %in% c("ageSpecRate", "ageSpecProj", "ageSpecBoth", "ageStdRate", "ageStdProj", "ageStdBoth")))
#     stop("type must be one of ageSpecRate, ageSpecProj, ageSpecBoth, ageStdRate, ageStdProj or ageStdBoth")
#   
#   if(type=="ageSpecBoth"){
#     projrate <- agespec.rate(APCList)
#     projcounts <- agespec.proj(APCList)
#     plab <- periodlabels(APCList)
#     
#     I <- nage(APCList)
#     J <- nperiod(APCList)
#     
#     y <- epi(APCList)
#     n <- pyrs(APCList)
#     
#     if(is.null(mfrow))
#       mfrow <- c(ceiling(I/4),4)
#     
#     if(is.null(coladd))
#       coladd=2
#     
#     par(mfrow=mfrow, ...)
#     for(u in 1:I){
#       if(showdata){
#         y_max <- max(y[,u]/n[,u]*scale, as.matrix(projrate[[u]])*scale, projcounts[[u]]/n[,u]*scale, na.rm=TRUE)
#       } else {
#         y_max <- max(projrate[[u]]*scale, projcounts[[u]]/n[,u]*scale)
#       } 
#       al <- agelabels(APCList)
#       matplot(plab, projrate[[u]]*scale,
#               main=paste("Age group", al[u]), type="n", lty=c(2,1,2), col=1,
#               xlab="Period", ylab=paste("Rate per ", format(scale, scientific=FALSE), sep=""), ylim=c(0, y_max), ...)
#       # polygon(c(plab, rev(plab)), c(projcounts[[u]][,1]/n[,u]*scale, rev(projcounts[[u]][,3]/n[,u]*scale)), col=coladd, border=NA)
#       polygon(c(plab, rev(plab)), c(projrate[[u]][,1]*scale, rev(projrate[[u]][,3]*scale)), col=coladd, border=NA)
#       matplot(plab, projcounts[[u]][,1:3]/n[,u]*scale, lty=c(2,1,2), col=1, add=T, type="l", ...)
#       lines(plab,  projrate[[u]][,2]*scale)
#       if(showdata){
#         if(is.na(sum(y))){
#           matplot(plab[1:(J-APCList@npred)], (y[,u]/n[,u]*scale)[1:(J-APCList@npred)], type="p", pch=20, add=T)
#         } else {
#           matplot(plab, y[,u]/n[,u]*scale, type="p", pch=20, add=T)
#         }
#       }
#       #abline(v=max(as.numeric(plab))-APCList@npred + .5, col="grey55", lwd=1.2)
#       numplab <- as.numeric(plab) 
#       abline(v=max(numplab)-(APCList@npred*(numplab[2]-numplab[1])) + .5, col="grey55", lwd=1.2)
#     }
#   }
#   
#   if(type=="ageSpecRate"){
#     projrate <- agespec.rate(APCList)
#     plab <- periodlabels(APCList)
#     
#     I <- nage(APCList)
#     J <- nperiod(APCList)
#     
#     y <- epi(APCList)
#     n <- pyrs(APCList)
#     
#     if(is.null(mfrow))
#       mfrow <- c(ceiling(I/4),4)
#     
#     par(mfrow=mfrow, ...)
#     for(u in 1:I){
#       if(showdata){
#         y_max <- max(y[,u]/n[,u]*scale, as.matrix(projrate[[u]][,1:3])*scale, na.rm=TRUE)
#       } else {
#         y_max <- max(projrate[[u]][,1:3]*scale)
#       } 
#       al <- agelabels(APCList)
#       matplot(plab, projrate[[u]][,1:3]*scale,
#               main=paste("Age group", al[u]), type="l", lty=c(2,1,2), col=1,
#               xlab="Period", ylab=paste("Rate per ",format(scale, scientific=FALSE), sep=""), ylim=c(0, y_max), ...)
#       if(showdata){
#         if(is.na(sum(y))){
#           matplot(plab[1:(J-APCList@npred)], (y[,u]/n[,u]*scale)[1:(J-APCList@npred)], type="p", pch=20, add=T)
#         } else {
#           matplot(plab, y[,u]/n[,u]*scale, type="p", pch=20, add=T)
#         }
#       }
#       #abline(v=max(as.numeric(plab))-APCList@npred + .5, col="grey55", lwd=1.2)
#       numplab <- as.numeric(plab) 
#       abline(v=max(numplab)-(APCList@npred*(numplab[2]-numplab[1])) + .5, col="grey55", lwd=1.2)
#       
#     }
#   }
#   if(type=="ageSpecProj"){
#     projcounts <- agespec.proj(APCList)
#     plab <- periodlabels(APCList)
#     
#     I <- nage(APCList)
#     J <- nperiod(APCList)
#     
#     y <- epi(APCList)
#     n <- pyrs(APCList)
#     
#     if(is.null(mfrow))
#       mfrow <- c(ceiling(I/4),4)
#     
#     par(mfrow=mfrow, ...)
#     for(u in 1:I){
#       if(showdata){
#         y_max <- max(y[,u]/n[,u]*scale,projcounts[[u]]/n[,u]*scale, na.rm=TRUE)
#       } else {
#         y_max <- max(projcounts[[u]]/n[,u]*scale)
#       } 
#       al <- agelabels(APCList)
#       matplot(plab, projcounts[[u]][,1:3]/n[,u]*scale,
#               main=paste("Age group", al[u]), type="l", lty=c(2,1,2), col=1,
#               xlab="Period", ylab=paste("Rate per ", format(scale, scientific=FALSE), sep=""), ylim=c(0, y_max), ...)
#       if(showdata){
#         if(is.na(sum(y))){
#           matplot(plab[1:(J-APCList@npred)], (y[,u]/n[,u]*scale)[1:(J-APCList@npred)], type="p", pch=20, add=T)
#         } else {
#           matplot(plab, y[,u]/n[,u]*scale, type="p", pch=20, add=T)
#         }
#       }
#       
#       #abline(v=max(as.numeric(plab))-APCList@npred + .5, col="grey55", lwd=1.2)
#       numplab <- as.numeric(plab) 
#       abline(v=max(numplab)-(APCList@npred*(numplab[2]-numplab[1])) + .5, col="grey55", lwd=1.2)
#       
#     }
#   } 
#   if(type=="ageStdBoth"){
#     J <- nperiod(APCList)
#     plab <- periodlabels(APCList)
#     my.wm <- matrix(rep(stdweight(APCList), each=J), byrow=F, nrow=J)
#     agg.n <- rowSums(pyrs(APCList))
#     
#     if(is.null(coladd))
#       coladd=2
#     
#     matplot(plab, agestd.rate(APCList)[,1:3]*scale, type="l", lty=c(2,1,2), col=1,  
#             xlab="Period", ylab=paste("Agestd. rate per ", format(scale, scientific=FALSE), sep=""),  ...)
#     polygon(c(plab, rev(plab)), c(agestd.rate(APCList)[,1]*scale, rev(agestd.rate(APCList)[,3]*scale)), col=coladd, border=NA)
#     matplot(plab, agestd.proj(APCList)[,1:3]/agg.n*scale, type="l", lty=c(2,1,2), col=c(1,1,1),  
#             add=T, ...)
#     lines(plab,  agestd.rate(APCList)[,2]*scale)
#     if(showdata){
#       y <- epi(APCList)
#       n <- pyrs(APCList)
#       if(is.na(sum(y))){
#         matplot(plab[1:(J-APCList@npred)], rowSums(my.wm * y/n*scale)[1:(J-APCList@npred)], type="p", pch=20, add=T)
#       } else {
#         matplot(plab, rowSums(my.wm * y/n*scale), type="p", pch=20, add=T)
#       }
#       
#       #abline(v=max(as.numeric(plab))-APCList@npred + .5, col="grey55", lwd=1.2)
#       numplab <- as.numeric(plab) 
#       abline(v=max(numplab)-(APCList@npred*(numplab[2]-numplab[1])) + .5, col="grey55", lwd=1.2)
#       
#     }
#   }
#   
#   if(type=="ageStdRate"){
#     J <- nperiod(APCList)
#     plab <- periodlabels(APCList)
#     my.wm <- matrix(rep(stdweight(APCList), each=J), byrow=F, nrow=J)
#     
#     matplot(plab, agestd.rate(APCList)[,1:3]*scale, type="l", lty=c(2,1,2), col=1,  
#             xlab="Period", ylab=paste("Agestd. rate per ", format(scale, scientific=FALSE), sep=""), ...)
#     if(showdata){
#       y <- epi(APCList)
#       n <- pyrs(APCList)
#       if(is.na(sum(y))){
#         matplot(plab[1:(J-APCList@npred)], rowSums(my.wm * y/n*scale)[1:(J-APCList@npred)], type="p", pch=20, add=T)
#       } else {
#         matplot(plab, rowSums(my.wm * y/n*scale), type="p", pch=20, add=T)
#       }
#       #abline(v=max(as.numeric(plab))-APCList@npred + .5, col="grey55", lwd=1.2)
#       numplab <- as.numeric(plab) 
#       abline(v=max(numplab)-(APCList@npred*(numplab[2]-numplab[1])) + .5, col="grey55", lwd=1.2)
#       
#     }
#   }
#   if(type=="ageStdProj"){
#     J <- nperiod(APCList)
#     plab <- periodlabels(APCList)
#     my.wm <- matrix(rep(stdweight(APCList), each=J), byrow=F, nrow=J)
#     agg.n <- rowSums(pyrs(APCList))
#     matplot(plab, agestd.proj(APCList)[,1:3]/agg.n*scale, type="l", lty=c(2,1,2), col=1,  
#             xlab="Period", ylab=paste("Agestd. rate per ", format(scale, scientific=FALSE), sep=""), ...)
#     if(showdata){
#       y <- epi(APCList)
#       n <- pyrs(APCList)
#       if(is.na(sum(y))){
#         matplot(plab[1:(J-APCList@npred)], rowSums(my.wm * y/n*scale, na.rm=T)[1:(J-APCList@npred)], type="p", pch=20, add=T)
#       } else {
#         matplot(plab, rowSums(my.wm * y/n*scale, na.rm=T), type="p", pch=20, add=T)
#       }
#       #abline(v=max(as.numeric(plab))-APCList@npred + .5, col="grey55", lwd=1.2)
#       numplab <- as.numeric(plab) 
#       abline(v=max(numplab)-(APCList@npred*(numplab[2]-numplab[1])) + .5, col="grey55", lwd=1.2)
#       
#     }
#   }
# }
# 
