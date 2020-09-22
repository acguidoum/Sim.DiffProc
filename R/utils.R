## Thu Apr 30 10:50:58 2020
## Original file Copyright © 2020 A.C. Guidoum, K. Boukhetala
## This file is part of the R package Sim.DiffProc
## Department of Probabilities & Statistics
## Faculty of Mathematics
## University of Science and Technology Houari Boumediene
## BP 32 El-Alia, U.S.T.H.B, Algiers
## Algeria

## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## A copy of the GNU General Public License is available at
## http://www.r-project.org/Licenses/
## Unlimited use and distribution (see LICENCE).
###################################################################################################



###
### Computes the bound of the confidence

bconfint <- function(x, ...)  UseMethod("bconfint")

bconfint.default <- function(x,level = 0.95,...)
            {
   if (!is.numeric(x)) 
         stop("argument 'x' must be numeric")
   return(quantile(x,c(0.5*(1-level), 1-0.5*(1-level)),type=7,na.rm=TRUE,...))
         }
###

# add.bconfint <- function(x,...) UseMethod("add.bconfint")

# add.bconfint.default <- function(x,level = 0.95,...)
            # {
    # lines(time(x),bconfint(x,level)[,1],...)
    # lines(time(x),bconfint(x,level)[,2],...)
         # }

###
### Computes the moment

moment <- function(x, ...)  UseMethod("moment")

moment.default <- function(x, order = 1,center = TRUE,...)
            {
	if (!is.numeric(x)) 
	       stop("argument 'x' must be numeric")
    if (any(!is.numeric(order)  || (order - floor(order) > 0) || order < 1)) 
	       stop(" 'order' must be a positive integer ")
	x = x[!is.na(x)]
    if (center) x <- x - mean(x)
    return(sum(x^order)/length(x))
         }
###
###

Mode <- function(x, ...)  UseMethod("Mode")

Mode.default <- function(x,...)
            {
    if (!is.numeric(x)) 
	      stop("argument 'x' must be numeric")
	x = x[!is.na(x)]
    return(density.default(x,na.rm = TRUE)$x[which.max(density.default(x,na.rm = TRUE)$y)] )
         }
		 
###
###

Median <- function(x, ...)  UseMethod("Median")

Median.default <- function(x,...)
            {
    if (!is.numeric(x)) 
	     stop("argument 'x' must be numeric")
    return(median(x, na.rm = TRUE))
         }		 

###
### Measure of Relative Variability

cv <- function(x, ...)  UseMethod("cv")

cv.default <- function(x,...)
            {
		if (!is.numeric(x)) 
		     stop("argument 'x' must be numeric")
        return(sd(x,na.rm = TRUE)/mean(x,na.rm = TRUE,...))
         }

		 
###
### Computes the sample skewness

skewness <- function(x,...) UseMethod("skewness")

skewness.default <- function(x,...)
                   {
	 if (!is.numeric(x)) 
	       stop("argument 'x' must be numeric")
     return(mean((x-mean(x,na.rm = TRUE))^3,na.rm = TRUE)/sd(x,na.rm = TRUE)^3)
}

###
### Computes the sample kurtosis

kurtosis <- function(x,...) UseMethod("kurtosis")

kurtosis.default <- function(x,...)
                   {
	 if (!is.numeric(x)) 
	        stop("argument 'x' must be numeric")
     return(mean((x-mean(x,na.rm = TRUE))^4,na.rm = TRUE)/sd(x,na.rm = TRUE)^4)
}

###
### Plot 2D and 3D

plot2d   <- function(x,...) UseMethod("plot2d")
lines2d  <- function(x,...) UseMethod("lines2d")
points2d <- function(x,...) UseMethod("points2d")

plot2d.default <- function(x,...)
        {
    class(x) <- "plot2d"
    plot(x,...)
}

lines2d.default <- function(x,...)
        {
    class(x) <- "lines2d"
    lines(x,...)
}

points2d.default <- function(x,...)
        {
    class(x) <- "points2d"
    points(x,...)
}


plot3D   <- function(x, ...)  UseMethod("plot3D")

plot3DD <- function(X,Y,Z,display = c("persp","rgl"),col=NULL,lwd=NULL,pch=NULL,
                            type = NULL,cex=NULL,main=NULL,sub=NULL,xlab=NULL,ylab=NULL,
                            zlab=NULL,grid=NULL,angle=NULL,...)
                 {
    display <- match.arg(display)
     # if ((display == "rgl") && !(require(rgl)) ) 
                  # stop("The 'rgl' package is not available.")        
     # else if ((display == "persp") && !(require(scatterplot3d)) ) 
                  # stop("The 'scatterplot3d' package is not available.")

    if (is.null(lwd))  
	      lwd = 1
    if (is.null(col))  
	      col = 2
    if (is.null(pch))  
	      pch = 16
    if (is.null(cex))  
	      cex = 0.6
    if (is.null(type)) 
	      type = "l"
    if (is.null(main)) 
	      main = ""
    if (is.null(sub))  
	      sub = ""
    if (is.null(xlab)) 
	      xlab = expression(x)
    if (is.null(ylab)) 
	      ylab = expression(y)
    if (is.null(zlab)) 
	      zlab = expression(z) 
    if (is.null(grid)) 
	      grid = TRUE   
    if (is.null(angle))
	      angle = 140
		  
    if (display=="persp"){
         scatterplot3d::scatterplot3d(X,Y,Z,angle =angle,color=col,lwd=lwd,type=type,pch=pch,
             main = main, sub = sub,xlab = xlab, ylab = ylab, zlab = zlab,
             grid = grid,cex.symbols=cex,...)
    }else if (display=="rgl"){
         rgl::plot3d(X,Y,Z,col=col,lwd=lwd,type=type,main = main, sub = sub,
             xlab = xlab, ylab = ylab, zlab = zlab,size=cex,...)
         }
}

plot3D.default <- function(x,display = c("persp", "rgl"),...) plot3DD(x,display,...)

####
#### plot for calss snssde2d

.plot.snssde2d <- function(x,union = TRUE,legend=TRUE,pos=1,type="l",main=NULL,col=NULL,lty=NULL,lwd=NULL,cex=NULL,
                          las=NULL,text.col=NULL,...) 
              {
    class(x) <- "snssde2d"
    if (is.null(col))
	     col=c(1,2)
    if (is.null(lty))
	     lty=c(1,1)
    if (is.null(lwd))
	     lwd=c(1,1)
    if (is.null(cex)) 
	     cex=0.75
    if (is.null(las))
	     las=1
    if (is.null(main))
	     main=""
    if (is.null(text.col))
	     text.col=c(1,1)
    if (pos==1)
	     pos = "top"
    else if (pos==2)
	     pos = "topright"
    else if (pos==3)
	     pos = "topleft"
    else if (pos==4)
	     pos = "center"
    else if (pos==5)
	     pos = "right"
    else if (pos==6)
	     pos = "left"
    else if (pos==7)
	     pos = "bottom"
    else if (pos==8)
	     pos = "bottomright"
    else if (pos==9)
	     pos = "bottomleft"
		 
	if (type=="n"){legend=FALSE}
    if (union==TRUE){
    plot(ts.union(x$X,x$Y),plot.type="single",ylab="",type="n",las=las,main=main,...)
    if (x$M == 1){
    lines(time(x),x$X,col=col[1],lty=lty[1],lwd=lwd[1],type=type,...)
    lines(time(x),x$Y,col=col[2],lty=lty[2],lwd=lwd[2],type=type,...)}else{
    for (i in 1:x$M){
    lines(time(x),x$X[,i],col=col[1],lty=lty[1],lwd=lwd[1],type=type,...)
    lines(time(x),x$Y[,i],col=col[2],lty=lty[2],lwd=lwd[2],type=type,...)
         }
    }
    if (legend){	
    legend(pos,c(expression(X[t]),expression(Y[t])),inset = .01,col=col,lty=lty,lwd=lwd,cex=cex,text.col=text.col)}
    }else{
    plot(x$X,plot.type="single",ylab=expression(X[t]),col=col[1],lty=lty[1],lwd=lwd[1],las=las,main=main,type=type,...)
    dev.new()
    plot(x$Y,plot.type="single",ylab=expression(Y[t]),col=col[2],lty=lty[2],lwd=lwd[2],las=las,main=main,type=type,...)
    }
}

.plot2d.snssde2d <- function(x,type=NULL,xlab=NULL,ylab=NULL,...)
                 {
    class(x) <- "snssde2d"
    if (is.null(type))
	     type="l"
    if (is.null(ylab))
	     ylab=expression(y)
    if (is.null(xlab))
	     xlab=expression(x)
    if (x$M == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)
                  Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{
                  X = x$X
                  Y = x$Y} 
    plot2d(X[,1],Y[,1],type=type,ylab=ylab,xlab=xlab,...)
    for(i in 3:4) axis(i)
}

####
#### plot for calss snssde3d

.plot.snssde3d <- function(x,union = TRUE,legend=TRUE,pos=1,type="l",main=NULL,col=NULL,lty=NULL,lwd=NULL,cex=NULL,
                          las=NULL,text.col=NULL,...) 
              {
    class(x) <- "snssde3d"
    if (is.null(col)) 
	    col=c(1,2,3)
    if (is.null(lty))
	    lty=c(1,1,1)
    if (is.null(lwd))
	    lwd=c(1,1,1)
    if (is.null(cex))
	    cex=0.75
    if (is.null(las))
	     las=1
    if (is.null(main))
	     main=""
    if (is.null(text.col))
	     text.col=c(1,1,1)
    if (pos==1)
	     pos = "top"
    else if (pos==2)
	     pos = "topright"
    else if (pos==3)
	     pos = "topleft"
    else if (pos==4)
	     pos = "center"
    else if (pos==5)
	     pos = "right"
    else if (pos==6)
	     pos = "left"
    else if (pos==7)
	     pos = "bottom"
    else if (pos==8)
	     pos = "bottomright"
    else if (pos==9)
	     pos = "bottomleft"
    if (type=="n")
	     legend=FALSE
    if (union){
    plot(ts.union(x$X,x$Y,x$Z),plot.type="single",ylab="",type="n",las=las,main=main,...)
    if (x$M == 1){
    lines(time(x),x$X,col=col[1],lty=lty[1],lwd=lwd[1],type,...)
    lines(time(x),x$Y,col=col[2],lty=lty[2],lwd=lwd[2],type,...)
    lines(time(x),x$Z,col=col[3],lty=lty[3],lwd=lwd[3],type,...)}else{
    for (i in 1:x$M){
    lines(time(x),x$X[,i],col=col[1],lty=lty[1],lwd=lwd[1],type,...)
    lines(time(x),x$Y[,i],col=col[2],lty=lty[2],lwd=lwd[2],type,...)
    lines(time(x),x$Z[,i],col=col[3],lty=lty[3],lwd=lwd[3],type,...)
         }
    }
    if (legend){	
    legend(pos,c(expression(X[t]),expression(Y[t]),expression(Z[t])),inset = .01,col=col,lty=lty,lwd=lwd,cex=cex,text.col=text.col)}
    }else{
    plot(x$X,plot.type="single",ylab=expression(X[t]),col=col[1],lty=lty[1],lwd=lwd[1],las=las,main=main,type=type,...)
    dev.new()
    plot(x$Y,plot.type="single",ylab=expression(Y[t]),col=col[2],lty=lty[2],lwd=lwd[2],las=las,main=main,type=type,...)
    dev.new()
    plot(x$Z,plot.type="single",ylab=expression(Z[t]),col=col[3],lty=lty[3],lwd=lwd[3],las=las,main=main,type=type,...)
    }
}


####
#### plot for calss bridgesde2d 

.plot.bridgesde2d  <- function(x,union = TRUE,legend=TRUE,pos=1,main=NULL,col=NULL,lty=NULL,lwd=NULL,cex=NULL,
                          las=NULL,text.col=NULL,type="l",...) 
              {
    class(x) <- "bridgesde2d"
    if (is.null(col))
	     col=c(1,2)
    if (is.null(lty))
	     lty=c(1,1)
    if (is.null(lwd))
	     lwd=c(1,1)
    if (is.null(cex))
	     cex=0.75
    if (is.null(las))
	     las=1
    if (is.null(main))
	     main=""
    if (is.null(text.col))
	     text.col=c(1,1)
    if (pos==1)
	     pos = "top"
    else if (pos==2)
	     pos = "topright"
    else if (pos==3)
	     pos = "topleft"
    else if (pos==4)
	     pos = "center"
    else if (pos==5)
	     pos = "right"
    else if (pos==6)
	     pos = "left"
    else if (pos==7)
	     pos = "bottom"
    else if (pos==8)
	     pos = "bottomright"
    else if (pos==9)
	     pos = "bottomleft"
	if (type=="n")
	    legend=FALSE
    if (union){
    plot(ts.union(x$X,x$Y),plot.type="single",ylab="",type="n",las=las,main=main,...)
    if (length(which(!is.na(x$Cx))) == 1 ){lines(as.vector(time(x$X)),x$X,col=col[1],lty=lty[1],lwd=lwd[1],type,...)}else{
	for (i in 1:length(which(!is.na(x$Cx)))) {lines(as.vector(time(x$X)),x$X[,i],col=col[1],lty=lty[1],lwd=lwd[1],type,...)}}
    if (length(which(!is.na(x$Cy))) == 1 ){lines(as.vector(time(x$Y)),x$Y,col=col[2],lty=lty[2],lwd=lwd[2],type,...)}else{
	for (i in 1:length(which(!is.na(x$Cy)))) {lines(as.vector(time(x$Y)),x$Y[,i],col=col[2],lty=lty[2],lwd=lwd[2],type,...)}}
    if (legend){	
    legend(pos,c(expression(X[t]),expression(Y[t])),inset = .01,col=col,lty=lty,lwd=lwd,cex=cex,text.col=text.col)}
    }else{
    plot(x$X,plot.type="single",ylab=expression(X[t]),col=col[1],lty=lty[1],lwd=lwd[1],las=las,main=main,type=type,...)
    dev.new()
    plot(x$Y,plot.type="single",ylab=expression(Y[t]),col=col[2],lty=lty[2],lwd=lwd[2],las=las,main=main,type=type,...)
    }
}

.plot2d.bridgesde2d <- function(x,type=NULL,xlab=NULL,ylab=NULL,...)
                 {
    class(x) <- "bridgesde2d"
    if (is.null(type))
	      type="l"
    if (is.null(ylab))
	      ylab=expression(Y[t])
    if (is.null(xlab))
	      xlab=expression(X[t])
    if (length(which(!is.na(x$Cx))) == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}  
    if (length(which(!is.na(x$Cy))) == 1){Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{Y = x$Y}	
    plot2d(X[,1],Y[,1],type=type,ylab=ylab,xlab=xlab,...)
    for(i in 3:4) axis(i)
}

####
#### plot for calss bridgesde3d
 
.plot.bridgesde3d <- function(x,union = TRUE,legend=TRUE,pos=1,type="l",main=NULL,col=NULL,lty=NULL,lwd=NULL,cex=NULL,
                          las=NULL,text.col=NULL,...) 
              {
    class(x) <- "snssde3d"
    if (is.null(col))
	      col=c(1,2,3)
    if (is.null(lty))
	      lty=c(1,1,1)
    if (is.null(lwd))
	      lwd=c(1,1,1)
    if (is.null(cex)) 
	      cex=0.75
    if (is.null(las))
	      las=1
    if (is.null(main))
	      main=""
    if (is.null(text.col))
	      text.col=c(1,1,1)
    if (pos==1)
	     pos = "top"
    else if (pos==2)
	     pos = "topright"
    else if (pos==3)
	     pos = "topleft"
    else if (pos==4)
	     pos = "center"
    else if (pos==5)
	     pos = "right"
    else if (pos==6)
	     pos = "left"
    else if (pos==7)
	     pos = "bottom"
    else if (pos==8)
	     pos = "bottomright"
    else if (pos==9)
	     pos = "bottomleft"
    if (type=="n")
	    legend=FALSE
    if (union){
    plot(ts.union(x$X,x$Y,x$Z),plot.type="single",ylab="",type="n",las=las,main=main,...)
    if (length(which(!is.na(x$Cx))) == 1 ){lines(as.vector(time(x$X)),x$X,col=col[1],lty=lty[1],lwd=lwd[1],type,...)}else{
	for (i in 1:length(which(!is.na(x$Cx)))) {lines(as.vector(time(x$X)),x$X[,i],col=col[1],lty=lty[1],lwd=lwd[1],type,...)}}
    if (length(which(!is.na(x$Cy))) == 1 ){lines(as.vector(time(x$Y)),x$Y,col=col[2],lty=lty[2],lwd=lwd[2],type,...)}else{
	for (i in 1:length(which(!is.na(x$Cy)))) {lines(as.vector(time(x$Y)),x$Y[,i],col=col[2],lty=lty[2],lwd=lwd[2],type,...)}}
    if (length(which(!is.na(x$Cz))) == 1 ){lines(as.vector(time(x$Z)),x$Z,col=col[3],lty=lty[3],lwd=lwd[3],type,...)}else{
	for (i in 1:length(which(!is.na(x$Cz)))) {lines(as.vector(time(x$Z)),x$Z[,i],col=col[3],lty=lty[3],lwd=lwd[3],type,...)}}
    if (legend){	
    legend(pos,c(expression(X[t]),expression(Y[t]),expression(Z[t])),inset = .01,col=col,lty=lty,lwd=lwd,cex=cex,text.col=text.col)}
    }else{
    plot(x$X,plot.type="single",ylab=expression(X[t]),col=col[1],lty=lty[1],lwd=lwd[1],las=las,main=main,type=type,...)
    dev.new()
    plot(x$Y,plot.type="single",ylab=expression(Y[t]),col=col[2],lty=lty[2],lwd=lwd[2],las=las,main=main,type=type,...)
    dev.new()
    plot(x$Z,plot.type="single",ylab=expression(Z[t]),col=col[3],lty=lty[3],lwd=lwd[3],las=las,main=main,type=type,...)
    }
}


####
#### plot for calss fptsde1d

.plot.dfptsde1d <- function(x,hist=FALSE,dens = NULL,legend=TRUE,pos=2,main=NULL,xlab=NULL,ylab=NULL,col=NULL,border=NULL,
                          lwd=NULL,cex=NULL,las=NULL,lty=NULL,xlim=NULL,ylim=NULL,add=NULL,...) 
               {
    class(x) <- "dfptsde1d"
    if (is.null(col)){col= c(rgb(255,0,0,75,maxColorValue=255),rgb(0,0,255,75,maxColorValue=255))}
    if (is.null(border)){border = c(rgb(255,0,0,130,maxColorValue=255),rgb(0,0,255,130,maxColorValue=255))}
    if (is.null(lty)){lty=1}
    if (is.null(lwd)){lwd=1}
    if (is.null(cex)){cex=1}
    if (is.null(las)){las=1}
    if (is.null(main)){main=""}
	if (is.null(add)){add=FALSE}
    if (is.null(xlab)){xlab="Time"}
	if (is.null(ylab)){ylab="Density"}
    if (pos==1)
	     pos = "top"
    else if (pos==2)
	     pos = "topright"
    else if (pos==3)
	     pos = "topleft"
    else if (pos==4)
	     pos = "center"
    else if (pos==5)
	     pos = "right"
    else if (pos==6)
	     pos = "left"
    else if (pos==7)
	     pos = "bottom"
    else if (pos==8)
	     pos = "bottomright"
    else if (pos==9)
	     pos = "bottomleft"
	
	if(hist==FALSE){
    if (is.null(dens)){
	    if (is.null(xlim)) {xlim = c(min(c(x$res$x),na.rm=TRUE) , max(c(x$res$x),na.rm=TRUE))}
        if (is.null(ylim)) {ylim = c(min(c(x$res$y),na.rm=TRUE) , max(c(x$res$y),na.rm=TRUE))}
	if (add==FALSE){
    plot.default(x$res , type = "n" , ylab = ylab , xlab = xlab,main=main,las=las,xlim = xlim , ylim = ylim, ...)
    polygon(x$res$x,x$res$y, col = col[1] , border = border[1],lty=lty,lwd=lwd,...)}else{
	polygon(x$res$x,x$res$y, col = col[2] , border = border[2],lty=lty,lwd=lwd,...)}}else{
    if (is.null(xlim)) {xlim = c(min(c(x$res$x, dens(x$res$y)),na.rm=TRUE) , max(c(x$res$x, dens(x$res$y)),na.rm=TRUE))}
    if (is.null(ylim)) {ylim = c(min(c(x$res$y, dens(x$res$x)),na.rm=TRUE) , max(c(x$res$y, dens(x$res$x)),na.rm=TRUE))}
    plot.default(x$res , type = "n" , ylab = ylab , xlab = xlab,main=main,las=las,xlim = xlim , ylim = ylim, ...)
    polygon(x$res$x,x$res$y, col = col[1] , border = border[1],lty=lty,lwd=lwd,...)
    curve(dens,col=1,lwd=2,add=TRUE,n = 1001) 
    if (legend){
    legend(pos,legend=c("Estimated density","True density"), inset = .01,lty = c(lty, 1),pch=c(NA,NA),
           lwd=c(2,2),col=c(col[1],1),cex=cex)    } 
    }}else{
	if (is.null(dens)){
		#if (is.null(xlim)) {xlim = c(min(c(x$res$x),na.rm=TRUE) , max(c(x$res$x),na.rm=TRUE))}
        #if (is.null(ylim)) {ylim = c(min(c(x$res$y),na.rm=TRUE) , max(c(x$res$y),na.rm=TRUE))}
	MASS::truehist(x$ech,xlab = xlab,ylab=ylab,main=main,las=las,col=col[1],...);box()}else{
	if (is.null(xlim)) {xlim = c(min(c(x$ech, dens(x$ech)),na.rm=TRUE) , max(c(x$ech, dens(x$ech)),na.rm=TRUE))}
    if (is.null(ylim)) {ylim = c(min(c(x$res$y, dens(x$res$x)),na.rm=TRUE) , max(c(x$res$y, dens(x$res$x)),na.rm=TRUE))}
	MASS::truehist(x$ech,xlab = xlab,ylab=ylab,main=main,las=las,col=col[1],xlim = xlim , ylim = ylim,...);box()
	curve(dens,col=1,lwd=2,add=TRUE,n = 1001)
	if (legend){
    legend(pos,legend=c("Distribution histogram","True density"), inset = .01,lty = c(NA, 1),pch=c(15,NA),
           lwd=c(2,2),col=c(col[1],1),cex=cex)    } 
	}
	}
}

####
#### plot for calss fptsde2d

.plot.dfptsde2d <- function(x,display=c("persp","rgl","image","contour"),hist=FALSE,drawpoints=FALSE,legend=TRUE,pos=2,main=NULL,xlab=NULL,ylab=NULL,zlab=NULL,col.pt=NULL,color.palette=NULL,col=NULL,border=NULL,
                          lwd=NULL,cex=NULL,pch=NULL,las=NULL,lty=NULL,xlim=NULL,ylim=NULL,expand = NULL,phi=NULL,theta=NULL,
                          ltheta=NULL,ticktype=NULL,shade=NULL,add=NULL,addline=FALSE,col2d=NULL,...) 
{
    class(x) <- "dfptsde2d"
    if (is.null(col)){col= c(rgb(255,0,0,75,maxColorValue=255),rgb(0,0,255,75,maxColorValue=255))}
    if (is.null(border)){border = c(rgb(255,0,0,130,maxColorValue=255),rgb(0,0,255,130,maxColorValue=255))}
	if (is.null(color.palette)){color.palette=colorRampPalette(c('white','green','blue','orange','red'))}
    if (is.null(col2d)){col2d = "lightblue"}
    if (is.null(lty)){lty=1}
    if (is.null(lwd)){lwd=1}
    if (is.null(cex)){cex=1}
    if (is.null(las)){las=1}
    if (is.null(pch)){pch=1}
    if (is.null(main)){main=""}
    if (is.null(col.pt)){col.pt="blue"}
    if (is.null(xlab)){xlab="x"}
    if (is.null(ylab)){ylab="y"}
    if (is.null(zlab)){zlab="Density function"}
    if (is.null(expand)){expand = .5}
    if (is.null(theta)) {theta = 30}
    if (is.null(ltheta)) {ltheta = 120}
    if (is.null(phi)) {phi = 30}
    if (is.null(shade)) {shade=0.75}
    if (is.null(ticktype)) {ticktype="detailed"}
    if (is.null(add)){add=FALSE}
    if (pos==1)
	     pos = "top"
    else if (pos==2)
	     pos = "topright"
    else if (pos==3)
	     pos = "topleft"
    else if (pos==4)
	     pos = "center"
    else if (pos==5)
	     pos = "right"
    else if (pos==6)
	     pos = "left"
    else if (pos==7)
	     pos = "bottom"
    else if (pos==8)
	     pos = "bottomright"
    else if (pos==9)
	     pos = "bottomleft"
    display <- match.arg(display)
    if (x$pdf =='Marginal'){
         if (hist==FALSE){ 
		if (missing(xlim)) {xlim = c(min(c(x$resx$x, x$resy$x)) , max(c(x$resx$x, x$resy$x)))}
		if (missing(ylim)) {ylim = c(min(c(x$resx$y, x$resy$y)) , max(c(x$resx$y, x$resy$y)))}
		plot(0 , type = "n" , ylab = "Density" , xlab = "Time" , las = 1 , xlim = xlim , ylim = ylim,main=main, ...)
		polygon(x$resx$x , x$resx$y, col = col[1] , border = border[1])
		polygon(x$resy$x , x$resy$y, col = col[2] , border = border[2])
         if (legend){	
         legend(pos,c(expression(hat(f)(x)),expression(hat(f)(y))),inset = .01,col=col,lty=lty,lwd=lwd,cex=cex)
         }}else{
	      MASS::truehist(x$fpt$x,xlab = "Time",main=main,las=las,col=col[1],...);box()
          dev.new()
          MASS::truehist(x$fpt$y,xlab = "Time",main=main,las=las,col=col[2],...);box()
    }}else if (x$pdf =='Joint'){
    #col2 <- heat.colors(length(x$res$z))[rank(x$res$z)]
    if (display=="persp"){
      persp(x$res$x,x$res$y,x$res$z,col=col2d,xlab=xlab,ylab=ylab,zlab=zlab,main=main,expand = expand,theta=theta,phi=phi,ltheta=ltheta,
      shade=shade,ticktype=ticktype,...)
    }else if (display=="rgl"){
      rgl::persp3d(x$res$x,x$res$y,x$res$z,col = col2d,xlab=xlab,ylab=ylab,zlab=zlab,main=main,...)
	  if(addline==TRUE) rgl::persp3d(x$res$x,x$res$y,x$res$z, front = "lines", back = "lines", lit = FALSE, add = TRUE)
    }else if (display=="image"){  
      image(x$res$x,x$res$y,x$res$z, col = col,xlab=xlab,ylab=ylab,main=main,las=las,add=add,...)
      contour(x$res$x,x$res$y,x$res$z, add = TRUE,...);box()
    if (drawpoints) {points(x$fpt$x, x$fpt$y, col=col.pt, cex=cex, pch=pch)}
    }else if (display=="contour"){
      filled.contour(x$res$x,x$res$y,x$res$z,plot.title = title(main = main,xlab = xlab, ylab = ylab),color.palette=color.palette,...)
    }
    }
}

####
#### plot for calss fptsde3d

.plot.dfptsde3d <- function(x,display="rgl",hist=FALSE,legend=TRUE,pos=2,main=NULL,xlab=NULL,ylab=NULL,zlab=NULL,col=NULL,border=NULL,
                          lwd=NULL,cex=NULL,las=NULL,lty=NULL,xlim=NULL,ylim=NULL,color.palette=NULL,...) 
{
    class(x) <- "dfptsde3d"
    if (is.null(col)){col= c(rgb(255,0,0,75,maxColorValue=255),rgb(0,0,255,75,maxColorValue=255),rgb(0,255,0,75,maxColorValue=255))}
    if (is.null(border)){border = c(rgb(255,0,0,130,maxColorValue=255),rgb(0,0,255,130,maxColorValue=255),rgb(0,255,0,130,maxColorValue=255))}
	if (is.null(color.palette)){color.palette=colorRampPalette(c('white','green','blue','orange','red'))}
    if (is.null(lty)){lty=1}
    if (is.null(lwd)){lwd=1}
    if (is.null(cex)){cex=1}
    if (is.null(las)){las=1}
    if (is.null(main)){main=""}
    if (is.null(xlab)){xlab="x"}
    if (is.null(ylab)){ylab="y"}
    if (is.null(zlab)){zlab="z"}
    if (pos==1)
	     pos = "top"
    else if (pos==2)
	     pos = "topright"
    else if (pos==3)
	     pos = "topleft"
    else if (pos==4)
	     pos = "center"
    else if (pos==5)
	     pos = "right"
    else if (pos==6)
	     pos = "left"
    else if (pos==7)
	     pos = "bottom"
    else if (pos==8)
	     pos = "bottomright"
    else if (pos==9)
	     pos = "bottomleft"
    #display <- match.arg(display)
    if (x$pdf =='Marginal'){
        if (hist==FALSE){ 
		if (missing(xlim)) {xlim = c(min(c(x$resx$x, x$resy$x,x$resz$x)) , max(c(x$resx$x, x$resy$x, x$resz$x)))}
		if (missing(ylim)) {ylim = c(min(c(x$resx$y, x$resy$y,x$resz$y)) , max(c(x$resx$y, x$resy$y, x$resz$y)))}
		plot(0 , type = "n" , ylab = "Density" , xlab = "Time" , las = 1 , xlim = xlim , ylim = ylim,main=main, ...)
		polygon(x$resx$x , x$resx$y, col = col[1] , border = border[1])
		polygon(x$resy$x , x$resy$y, col = col[2] , border = border[2])
		polygon(x$resz$x , x$resz$y, col = col[3] , border = border[3])
         if (legend){	
         legend(pos,c(expression(hat(f)(x)),expression(hat(f)(y)),expression(hat(f)(z))),inset = .01,col=col,lty=lty,lwd=lwd,cex=cex)
         }}else{
	      MASS::truehist(x$fpt$x,xlab = "Time",main=main,las=las,col=col[1],...);box()
          dev.new()
          MASS::truehist(x$fpt$y,xlab = "Time",main=main,las=las,col=col[2],...);box()
          dev.new()
          MASS::truehist(x$fpt$z,xlab = "Time",main=main,las=las,col=col[3],...);box()
    }}else if (x$pdf =='Joint'){
      sm::sm.density(x$fpt,display="rgl",h=x$res$h,...)
    }
}

####
#### plot for calss dsde1d

.plot.dsde1d <- function(x,hist=FALSE,dens = NULL,legend=TRUE,pos=2,main=NULL,xlab=NULL,ylab=NULL,col=NULL,border=NULL,
                          lwd=NULL,cex=NULL,las=NULL,lty=NULL,xlim=NULL,ylim=NULL,add=NULL,...) 
               {
    class(x) <- "dsde1d"
    if (is.null(col)){col= c(rgb(255,0,0,75,maxColorValue=255),rgb(0,0,255,75,maxColorValue=255))}
    if (is.null(border)){border = c(rgb(255,0,0,130,maxColorValue=255),rgb(0,0,255,130,maxColorValue=255))}
    # if (is.null(col)){col= rgb(255,0,0,75,maxColorValue=255)}
    # if (is.null(border)){border = rgb(255,0,0,130,maxColorValue=255)}
    if (is.null(lty)){lty=1}
    if (is.null(lwd)){lwd=1}
    if (is.null(cex)){cex=1}
    if (is.null(las)){las=1}
    if (is.null(add)){add=FALSE}
    if (is.null(main)){main=""}
    if (is.null(xlab)){xlab="x"}
	if (is.null(ylab)){ylab="Density"}
    if (pos==1)
	     pos = "top"
    else if (pos==2)
	     pos = "topright"
    else if (pos==3)
	     pos = "topleft"
    else if (pos==4)
	     pos = "center"
    else if (pos==5)
	     pos = "right"
    else if (pos==6)
	     pos = "left"
    else if (pos==7)
	     pos = "bottom"
    else if (pos==8)
	     pos = "bottomright"
    else if (pos==9)
	     pos = "bottomleft"
	
	if(hist==FALSE){
    if (is.null(dens)){
	    if (is.null(xlim)) {xlim = c(min(c(x$res$x),na.rm=TRUE) , max(c(x$res$x),na.rm=TRUE))}
        if (is.null(ylim)) {ylim = c(min(c(x$res$y),na.rm=TRUE) , max(c(x$res$y),na.rm=TRUE))}
	if (add==FALSE){
    plot.default(x$res , type = "n" , ylab = ylab , xlab = xlab,main=main,las=las,xlim = xlim , ylim = ylim, ...)
    polygon(x$res$x,x$res$y, col = col[1] , border = border[1],lty=lty,lwd=lwd,...)}else{
	polygon(x$res$x,x$res$y, col = col[2] , border = border[2],lty=lty,lwd=lwd,...)}}else{
    if (is.null(xlim)) {xlim = c(min(c(x$res$x, dens(x$res$y)),na.rm=TRUE) , max(c(x$res$x, dens(x$res$y)),na.rm=TRUE))}
    if (is.null(ylim)) {ylim = c(min(c(x$res$y, dens(x$res$x)),na.rm=TRUE) , max(c(x$res$y, dens(x$res$x)),na.rm=TRUE))}
    plot.default(x$res , type = "n" , ylab = ylab , xlab = xlab,main=main,las=las,xlim = xlim , ylim = ylim, ...)
    polygon(x$res$x,x$res$y, col = col[1] , border = border[1],lty=lty,lwd=lwd,...)
    curve(dens,col=1,lwd=2,add=TRUE,n = 1001) 
    if (legend){
    legend(pos,legend=c("Estimated density","True density"), inset = .01,lty = c(lty, 1),pch=c(NA,NA),
           lwd=c(2,2),col=c(col[1],1),cex=cex)    } 
    }}else{
	if (is.null(dens)){
		#if (is.null(xlim)) {xlim = c(min(c(x$res$x),na.rm=TRUE) , max(c(x$res$x),na.rm=TRUE))}
        #if (is.null(ylim)) {ylim = c(min(c(x$res$y),na.rm=TRUE) , max(c(x$res$y),na.rm=TRUE))}
	MASS::truehist(x$ech,xlab = xlab,ylab=ylab,main=main,las=las,col=col[1],...);box()}else{
	if (is.null(xlim)) {xlim = c(min(c(x$ech, dens(x$ech)),na.rm=TRUE) , max(c(x$ech, dens(x$ech)),na.rm=TRUE))}
    if (is.null(ylim)) {ylim = c(min(c(x$res$y, dens(x$res$x)),na.rm=TRUE) , max(c(x$res$y, dens(x$res$x)),na.rm=TRUE))}
	MASS::truehist(x$ech,xlab = xlab,ylab=ylab,main=main,las=las,col=col[1],xlim = xlim , ylim = ylim,...);box()
	curve(dens,col=1,lwd=2,add=TRUE,n = 1001)
	if (legend){
    legend(pos,legend=c("Distribution histogram","True density"), inset = .01,lty = c(NA, 1),pch=c(15,NA),
           lwd=c(2,2),col=c(col[1],1),cex=cex)    } 
	}
	}
}


####
#### plot for calss dsde2d

.plot.dsde2d <- function(x,display=c("persp","rgl","image","contour"),hist=FALSE,drawpoints=FALSE,legend=TRUE,pos=2,main=NULL,xlab=NULL,ylab=NULL,zlab=NULL,col.pt=NULL,color.palette=NULL,col=NULL,border=NULL,
                          lwd=NULL,cex=NULL,pch=NULL,las=NULL,lty=NULL,xlim=NULL,ylim=NULL,expand = NULL,phi=NULL,theta=NULL,
                          ltheta=NULL,ticktype=NULL,shade=NULL,add=NULL,cex.main=NULL,addline=FALSE,col2d=NULL,...) 
{
    class(x) <- "dsde2d"
    if (is.null(col)){col= c(rgb(255,0,0,75,maxColorValue=255),rgb(0,0,255,75,maxColorValue=255))}
    if (is.null(border)){border = c(rgb(255,0,0,130,maxColorValue=255),rgb(0,0,255,130,maxColorValue=255))}
	if (is.null(color.palette)){color.palette=colorRampPalette(c('white','green','blue','orange','red'))}
    if (is.null(col2d)){col2d = "lightblue"}	
    if (is.null(lty)){lty=1}
    if (is.null(lwd)){lwd=1}
    if (is.null(cex)){cex=1}
    if (is.null(cex.main)){cex.main=1}
    if (is.null(las)){las=1}
    if (is.null(pch)){pch=1}
    if (is.null(main)){main=""}
    if (is.null(col.pt)){col.pt="blue"}
    if (is.null(xlab)){xlab="x"}
    if (is.null(ylab)){ylab="y"}
    if (is.null(zlab)){zlab="Density function"}
    if (is.null(expand)){expand = .5}
    if (is.null(theta)) {theta = 30}
    if (is.null(ltheta)) {ltheta = 120}
    if (is.null(phi)) {phi = 30}
    if (is.null(shade)) {shade=0.75}
    if (is.null(ticktype)) {ticktype="detailed"}
    if (is.null(add)){add=FALSE}
    if (pos==1)
	     pos = "top"
    else if (pos==2)
	     pos = "topright"
    else if (pos==3)
	     pos = "topleft"
    else if (pos==4)
	     pos = "center"
    else if (pos==5)
	     pos = "right"
    else if (pos==6)
	     pos = "left"
    else if (pos==7)
	     pos = "bottom"
    else if (pos==8)
	     pos = "bottomright"
    else if (pos==9)
	     pos = "bottomleft"
		 
    display <- match.arg(display)
    if (x$pdf =='Marginal'){
         if (hist==FALSE){ 
		if (missing(xlim)) {xlim = c(min(c(x$resx$x, x$resy$x)) , max(c(x$resx$x, x$resy$x)))}
		if (missing(ylim)) {ylim = c(min(c(x$resx$y, x$resy$y)) , max(c(x$resx$y, x$resy$y)))}
		plot(0 , type = "n" , ylab = "Density" , xlab = "State variables" , las = 1 , xlim = xlim , ylim = ylim,main=main, ...)
		polygon(x$resx$x , x$resx$y, col = col[1] , border = border[1])
		polygon(x$resy$x , x$resy$y, col = col[2] , border = border[2])
         if (legend){	
         legend(pos,c(expression(hat(f)(x)),expression(hat(f)(y))),inset = .01,col=col,lty=lty,lwd=lwd,cex=cex)
         }}else{
	     MASS::truehist(x$ech$x,xlab = "x",main=main,las=las,col=col[1],...);box()
          dev.new()
          MASS::truehist(x$ech$y,xlab = "y",main=main,las=las,col=col[2],...);box()
    }}else if (x$pdf =='Joint'){
    #col2 <- heat.colors(length(x$res$z))[rank(x$res$z)]
    if (display=="persp"){
      persp(x$res$x,x$res$y,x$res$z,col=col2d,xlab=xlab,ylab=ylab,zlab=zlab,main=main,expand = expand,theta=theta,phi=phi,ltheta=ltheta,
      shade=shade,ticktype=ticktype,...)
    }else if (display=="rgl"){
      rgl::persp3d(x$res$x,x$res$y,x$res$z,col = col2d,xlab=xlab,ylab=ylab,zlab=zlab,main=main,...)
	  if(addline==TRUE) rgl::persp3d(x$res$x,x$res$y,x$res$z, front = "lines", back = "lines", lit = FALSE, add = TRUE)
    }else if (display=="image"){  
      image(x$res$x,x$res$y,x$res$z, col = col2d,xlab=xlab,ylab=ylab,main=main,las=las,add=add,...)
      contour(x$res$x,x$res$y,x$res$z, add = TRUE,...);box()
    if (drawpoints) {points(x$ech$x, x$ech$y, col=col.pt, cex=cex, pch=pch)}
    }else if (display=="contour"){
      filled.contour(x$res$x,x$res$y,x$res$z,plot.title = title(main = main,xlab = xlab, ylab = ylab,cex.main=cex.main),color.palette=color.palette,...)
    }
    }
}

####
#### plot for calss dsde3d
# .plotks.3d <- function(fhat,display=c("rgl","persp"), cont=c(25,50,75), abs.cont, approx.cont=TRUE, colors, col.fun, alphavec, size=3, col.pt="blue", add=FALSE, xlab=NULL, ylab=NULL, zlab=NULL, drawpoints=FALSE, alpha=1, box=TRUE, axes=TRUE, ...)

# {
  # compute contours
  # if (missing(abs.cont))
  # {
    # if (!is.null(fhat$cont))
      # {
        # cont.ind <- rep(FALSE, length(fhat$cont))
          # for (j in 1:length(cont))
            # cont.ind[which(cont[j] == 100-as.numeric(unlist(strsplit(names(fhat$cont),"%"))))] <- TRUE
          
        # if (all(!cont.ind))
          # hts <- contourLevels(fhat, prob=(100-cont)/100, approx=approx.cont)
        # else
          # hts <- fhat$cont[cont.ind]
      # }
    # else
      # hts <- contourLevels(fhat, prob=(100-cont)/100, approx=approx.cont)
  # }  
  # else
    # hts <- abs.cont
  
  # nc <- length(hts)
  
  # if (missing(colors)) colors <- rev(heat.colors(nc))
  # if (!missing(col.fun)) colors <- col.fun(nc)
  # if (is.null(xlab)) xlab <- "x"
  # if (is.null(ylab)) ylab <- "y"
  # if (is.null(zlab)) zlab <- "z"
  # if (missing(alphavec))
  # {
    # if (is.null(fhat$deriv.order)) alphavec <- seq(0.1,0.5,length=nc)
    # else alphavec <- c(rev(seq(0.1,0.4,length=round(nc/2))), seq(0.1,0.4,length=round(nc/2)))
  # }
 
  # fhat.eval.mean <- sapply(fhat$eval.points, mean)
  # display <- match.arg(display)
  # if (display=="rgl"){
  # if (drawpoints)
    # rgl::plot3d(fhat$x[,1],fhat$x[,2],fhat$x[,3], size=size, col=col.pt, alpha=alpha, xlab=xlab, ylab=ylab, zlab=zlab, add=add, box=FALSE, axes=FALSE, ...)
  # else
    # rgl::plot3d(fhat$x[,1],fhat$x[,2],fhat$x[,3], size=0, col="transparent", alpha=0, xlab=xlab, ylab=ylab, zlab=zlab, add=add, box=FALSE, axes=FALSE, ...)  
  
  # for (i in 1:nc)
    # if (hts[nc-i+1] < max(fhat$estimate))
      # misc3d::contour3d(fhat$estimate, level=hts[nc-i+1], x=fhat$eval.points[[1]], y=fhat$eval.points[[2]], z=fhat$eval.points[[3]], add=TRUE, color=colors[i], alpha=alphavec[i], box=FALSE, axes=FALSE, ...)

  # if (axes) axes3d()
  # if (box) box3d()}else{
  #if (drawpoints)
  #scatterplot3d::scatterplot3d(fhat$x[,1],fhat$x[,2],fhat$x[,3], type="p",xlab="x",ylab="y",zlab="z",color=col.pt,axis=axes,box=box,xlim=c(min(fhat$eval.points[[1]],na.rm =TRUE),max(fhat$eval.points[[1]]),na.rm =TRUE),ylim=c(min(fhat$eval.points[[2]],na.rm =TRUE),max(fhat$eval.points[[2]]),na.rm =TRUE),zlim=c(min(fhat$eval.points[[3]],na.rm =TRUE),max(fhat$eval.points[[3]]),na.rm =TRUE),...)
  #else
  #scatterplot3d::scatterplot3d(fhat$x[,1],fhat$x[,2],fhat$x[,3], type="n",xlab="x",ylab="y",zlab="z",axis=axes,box=box,xlim=c(min(fhat$eval.points[[1]],na.rm =TRUE),max(fhat$eval.points[[1]],na.rm =TRUE)),ylim=c(min(fhat$eval.points[[2]],na.rm =TRUE),max(fhat$eval.points[[2]],na.rm =TRUE)),zlim=c(min(fhat$eval.points[[3]],na.rm =TRUE),max(fhat$eval.points[[3]],na.rm =TRUE)),...)
  # plot(0,type="n",axes=FALSE,xlab="",ylab="",xlim=c(min(fhat$eval.points[[1]],na.rm =TRUE),max(fhat$eval.points[[1]],na.rm =TRUE)),ylim=c(min(fhat$eval.points[[2]],na.rm =TRUE),max(fhat$eval.points[[2]],na.rm =TRUE)))
  # misc3d::contour3d(fhat$estimate, level=hts[nc-nc+1], x=fhat$eval.points[[1]], y=fhat$eval.points[[2]], z=fhat$eval.points[[3]], add=TRUE, color=colors[1], alpha=alphavec[nc], engine = "standard",perspective=0,...)
  # for (i in 1:(nc-1))
    # if (hts[nc-i+1] < max(fhat$estimate))
      # misc3d::contour3d(fhat$estimate, level=hts[nc-i+1], x=fhat$eval.points[[1]], y=fhat$eval.points[[2]], z=fhat$eval.points[[3]], add=TRUE, color=colors[i], alpha=alphavec[i], engine = "standard",perspective=0, ...)
# }
# }

.plot.dsde3d <- function(x,display="rgl",hist=FALSE,legend=TRUE,pos=2,main=NULL,xlab=NULL,ylab=NULL,zlab=NULL,col=NULL,border=NULL,
                          lwd=NULL,cex=NULL,las=NULL,lty=NULL,xlim=NA,ylim=NA,zlim=NA,...) 
{
    class(x) <- "dsde3d"
    if (is.null(col)){col= c(rgb(255,0,0,75,maxColorValue=255),rgb(0,0,255,75,maxColorValue=255),rgb(0,255,0,75,maxColorValue=255))}
    if (is.null(border)){border = c(rgb(255,0,0,130,maxColorValue=255),rgb(0,0,255,130,maxColorValue=255),rgb(0,255,0,130,maxColorValue=255))}
    if (is.null(lty)){lty=1}
    if (is.null(lwd)){lwd=1}
    if (is.null(cex)){cex=1}
    if (is.null(las)){las=1}
    if (is.null(main)){main=""}
    if (is.null(xlab)){xlab="x"}
    if (is.null(ylab)){ylab="y"}
    if (is.null(zlab)){zlab="z"}
    if (pos==1)
	     pos = "top"
    else if (pos==2)
	     pos = "topright"
    else if (pos==3)
	     pos = "topleft"
    else if (pos==4)
	     pos = "center"
    else if (pos==5)
	     pos = "right"
    else if (pos==6)
	     pos = "left"
    else if (pos==7)
	     pos = "bottom"
    else if (pos==8)
	     pos = "bottomright"
    else if (pos==9)
	     pos = "bottomleft"
		 
    if (x$pdf =="Marginal"){
        if (hist==FALSE){ 
		if (is.na(xlim)) {xlim = c(min(c(x$resx$x, x$resy$x,x$resz$x)) , max(c(x$resx$x, x$resy$x, x$resz$x)))}
		if (is.na(ylim)) {ylim = c(min(c(x$resx$y, x$resy$y,x$resz$y)) , max(c(x$resx$y, x$resy$y, x$resz$y)))}
		plot(0 , type = "n" , ylab = "Density" , xlab = "State variables" , las = 1 , xlim = xlim , ylim = ylim,main=main, ...)
		polygon(x$resx$x , x$resx$y, col = col[1] , border = border[1])
		polygon(x$resy$x , x$resy$y, col = col[2] , border = border[2])
		polygon(x$resz$x , x$resz$y, col = col[3] , border = border[3])
         if (legend){	
         legend(pos,c(expression(hat(f)(x)),expression(hat(f)(y)),expression(hat(f)(z))),inset = .01,col=col,lty=lty,lwd=lwd,cex=cex)
         }}else{
	     MASS::truehist(x$ech$x,xlab = "x",main=main,las=las,col=col[1],...);box()
          dev.new()
          MASS::truehist(x$ech$y,xlab = "y",main=main,las=las,col=col[2],...);box()
          dev.new()
          MASS::truehist(x$ech$z,xlab = "z",main=main,las=las,col=col[3],...);box()
    }}else if (x$pdf =="Joint"){
          sm::sm.density(x$ech,display="rgl",h=x$res$h,xlim=xlim,ylim=ylim,zlim=zlim,...)   
    }
}


####
#### Char2expression

.Char2exp <- function(expr)
         {
y <- gsub(pattern = 'theta[1]', replacement = 'theta1', x = expr, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[2]', replacement = 'theta2', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[3]', replacement = 'theta3', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[4]', replacement = 'theta4', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[5]', replacement = 'theta5', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[6]', replacement = 'theta6', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[7]', replacement = 'theta7', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[8]', replacement = 'theta8', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[9]', replacement = 'theta9', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[10]', replacement = 'theta10', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[11]', replacement = 'theta11', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[12]', replacement = 'theta12', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[13]', replacement = 'theta13', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[14]', replacement = 'theta14', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[15]', replacement = 'theta15', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[16]', replacement = 'theta16', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[17]', replacement = 'theta17', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[18]', replacement = 'theta18', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[19]', replacement = 'theta19', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[20]', replacement = 'theta20', x = y, ignore.case = F,fixed = T)
Y <- parse(text=y,srcfile=y)
return(Y)
}


.Exp2char <- function(expr)
         {
y <- gsub(pattern = 'theta1', replacement = 'theta[1]', x = expr, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta2', replacement = 'theta[2]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta3', replacement = 'theta[3]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta4', replacement = 'theta[4]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta5', replacement = 'theta[5]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta6', replacement = 'theta[6]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta7', replacement = 'theta[7]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta8', replacement = 'theta[8]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta9', replacement = 'theta[9]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta10', replacement = 'theta[10]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta11', replacement = 'theta[11]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta12', replacement = 'theta[12]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta13', replacement = 'theta[13]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta14', replacement = 'theta[14]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta15', replacement = 'theta[15]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta16', replacement = 'theta[16]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta17', replacement = 'theta[17]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta18', replacement = 'theta[18]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta19', replacement = 'theta[19]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta20', replacement = 'theta[20]', x = y, ignore.case = F,fixed = T)
Y <- parse(text=y,srcfile=y)
return(Y)
}

#####
##### Simulate a Multivariate Normal Distribution

# .rMnorm <- function(n , Mu, Sigma) 
         # {
    # d <- dim(Sigma)[1]
    # Egn <- eigen(Sigma)
    # Y <- matrix(rnorm(d * n), n)
    # X <- drop(Mu) + Egn$vectors %*% diag(sqrt(pmax(Egn$values, 0)), d) %*% t(Y)
    # return(t(X))
# }

#####
##### dc.sde1d

.dcEuler <- function(x, t, x0, t0, theta, drift, diffusion, log=FALSE)
              {
	 A   <- function(t,x,theta)  eval(drift)
     S   <- function(t,x,theta)  eval(diffusion)
     dt  <- t-t0
     lik <- dnorm(x,mean=x0+A(t0,x0,theta)*dt,sd=sqrt(dt)*S(t0,x0,theta),log=log)
     lik[is.infinite(lik)] <- NA
     lik
       }
	   
.dcOzaki <- function(x, t, x0, t0, theta, drift, diffusion, log=FALSE)
             {  
     A   <- function(t,x,theta)  eval(drift)
     S   <- function(t,x,theta)  eval(diffusion)
     Ax  <- function(t,x,theta)  eval(.Exp2char(as.expression(D(.Char2exp(drift),"x"))))
     dt  <- t-t0
     K   <- log(1+(A(t0,x0,theta)/(x0*Ax(t0,x0,theta)))*(exp(Ax(t0,x0,theta)*dt)-1))/dt
     E   <- x0 + (A(t0,x0,theta)/Ax(t0,x0,theta))*(exp(Ax(t0,x0,theta)*dt)-1)
     V   <- S(t0,x0,theta)^2 * ((exp(2*K*dt) -1)/(2*K))
     lik <- dnorm(x, mean=E, sd=sqrt(V),log=log)
     lik[is.infinite(lik)] <- NA
     lik      
}

.dcShoji <- function(x, t, x0, t0, theta, drift, diffusion, log=FALSE)
            {
     A   <- function(t,x,theta)  eval(drift)
     S   <- function(t,x,theta)  eval(diffusion)
	 At  <- function(t,x,theta)  eval(.Exp2char(as.expression(D(.Char2exp(drift),"t"))))
	 Ax  <- function(t,x,theta)  eval(.Exp2char(as.expression(D(.Char2exp(drift),"x"))))
	 Axx <- function(t,x,theta)  eval(.Exp2char(as.expression(D(D(.Char2exp(drift),"x"),"x"))))
     dt  <- t-t0
     E   <- x0 + A(t0,x0,theta)*(exp(Ax(t0,x0,theta)*dt)-1)/Ax(t0,x0,theta) + (0.5 * S(t0,x0,theta)^2 * Axx(t0,x0,theta) + At(t0,x0,theta))*(exp(Ax(t0,x0,theta)*dt)-1-Ax(t0,x0,theta)*dt)/Ax(t0,x0,theta)^2 
     V   <- S(t0,x0,theta)^2 * (exp(2*Ax(t0,x0,theta)*dt)-1)/(2*Ax(t0,x0,theta))
     lik <- dnorm(x, mean=E, sd=sqrt(V),log=log) 
	 lik[is.infinite(lik)] <- NA
     lik   
}

# .dcElerian <- function(x, t, x0, t0, theta, drift,diffusion, log=FALSE)
            # {
     # test <- .Exp2char(as.expression(D(.Char2exp(diffusion),"x")))[[1]]
     # if (test == 0) stop("The approximation is not valid, because 'deriv(diffusion,x) = 0'")
     # A   <- function(t,x,theta)  eval(drift)
     # S   <- function(t,x,theta)  eval(diffusion)
	 # Sx  <- function(t,x,theta)  eval(.Exp2char(as.expression(D(.Char2exp(diffusion),"x"))))	 
     # dt <- t-t0
     # E <- 0.5*S(t0, x0, theta)*Sx(t0, x0, theta)*dt
     # B <- -0.5* (S(t0, x0,theta)/Sx(t0, x0,theta)) + x0 + A(t0, x0, theta)*dt - E
     # C <- 1/((S(t0, x0,theta)^2) * dt)
     # z <- (x-B)/E
     # z[z<0] <- NA
     # lik <- ( exp(-0.5*(C+z)) * cosh(sqrt(C*z)) )/( sqrt(z) *abs(E)* sqrt(2*pi) )	 
     # if(log) lik <- log(lik)
     # lik[is.infinite(lik)] <- NA
     # lik
# }

.dcKessler <- function(x, t, x0, t0, theta, drift, diffusion, log=FALSE)
           {
     A   <- function(t,x,theta)  eval(drift)
	 Ax  <- function(t,x,theta)  eval(.Exp2char(as.expression(D(.Char2exp(drift),"x"))))
	 Axx <- function(t,x,theta)  eval(.Exp2char(as.expression(D(D(.Char2exp(drift),"x"),"x"))))	
     S   <- function(t,x,theta)  eval(diffusion)
	 Sx  <- function(t,x,theta)  eval(.Exp2char(as.expression(D(.Char2exp(diffusion),"x"))))
	 Sxx <- function(t,x,theta)  eval(.Exp2char(as.expression(D(D(.Char2exp(diffusion),"x"),"x"))))		   
     dt <- t-t0 
     E   <- x0 + A(t0,x0,theta)*dt + 0.5*(dt)^2 * (A(t0,x0,theta)*Ax(t0,x0,theta)+0.5*(S(t0,x0,theta))^2 * Axx(t0,x0,theta))
     V   <- x0^2 + (2*A(t0,x0,theta)*x0+(S(t0,x0,theta))^2)*dt + 0.5*(dt)^2 *(2*A(t0,x0,theta)*(Ax(t0,x0,theta)*x0+A(t0,x0,theta)+S(t0,x0,theta)*Sx(t0,x0,theta))+(S(t0,x0,theta))^2 *(Axx(t0,x0,theta)*x0+2*Ax(t0,x0,theta)+(Sx(t0,x0,theta))^2 + S(t0,x0,theta)*Sxx(t0,x0,theta))) - E^2
	 V[V < 0] <- NA
     lik <- dnorm(x, mean=E, sd=sqrt(V),log=log) 
	 lik[is.infinite(lik)] <- NA
     lik 
}
