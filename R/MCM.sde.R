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


#####
##### Monte-Carlo SDE's

MCM.sde <- function(model, ...)  UseMethod("MCM.sde")

MCM.sde.default <- function(model,statistic,R=100,time,exact=NULL,names=NULL,
                            level = 0.95, parallel = c("no", "multicore", "snow"),
                            ncpus = getOption("ncpus",1L), cl = NULL,...)
                 {
     if (!is.function(statistic)) 
	      stop("'statistic' must be function")
     if (!is.null(exact)){
        if (!is.list(exact)) 
		   stop("'exact' must be a named list")}
     if (!is.null(names)){
        if (!is.character(names)) 
		    stop("'names' must be an character")}
     if (any(!is.numeric(R)  || (R - floor(R) > 0) || R <= 0)) 
	        stop(" 'R' must be an integer ")
     if (R == 1) stop("Parallel Monte-Carlo Methods is for any R > 1")			
     if (missing(time)) {time = as.numeric(model$T)}
	 if (length(time) > 1 ) 
	        stop (" 'time' must be an integer 't0 < time <= T' ")
     if (any(model$T < time | model$t0 > time) )  
	        stop( " please use 't0 < time <= T'")
     if (any(level <= 0 | level >= 1) )  
	        stop( " please use '0 < level < 1'")
     if (class(model)=="bridgesde1d" | class(model)=="bridgesde2d" | class(model)=="bridgesde3d") 
	        stop("Not available for diffusion bridges")
     if (class(model)!="snssde1d" & class(model)!="snssde2d" & class(model)!="snssde3d") 
	        stop(" 'model' is not class of 'snssdekd', k=1,2,3. ")
     if (missing(parallel)) 
	     parallel <- getOption("parallel", "no")
     parallel <- match.arg(parallel)
     have_mc <- have_snow <- FALSE
     if (parallel != "no" && ncpus > 1L) {
         if (parallel == "multicore") 
            have_mc <- .Platform$OS.type != "windows"
         else if (parallel == "snow") 
            have_snow <- TRUE
         if (!have_mc && !have_snow) 
            ncpus <- 1L
        loadNamespace("parallel")
     }
     if (class(model)=="snssde1d"){
       statistic <- match.fun(statistic)
       fn1d      <- function(x,...) statistic(x, ...)
       rand      <- matrix(0,nrow=R,ncol=model$M)
       run1d     <- function(...) Sim.DiffProc::rsde1d(eval(model$call),at=time)
       if (ncpus > 1L && (have_mc || have_snow)){
         if (have_mc) {
		   base::RNGkind("L'Ecuyer-CMRG")
           rand <- do.call("rbind",parallel::mclapply(1:R, run1d,mc.cores = ncpus,mc.set.seed = TRUE))
         }else if (have_snow) {
            list(...)
            if (is.null(cl)){
                cl  <- parallel::makePSOCKcluster(rep("localhost", ncpus))
                parallel::clusterEvalQ(cl, library(Sim.DiffProc))
                parallel::clusterExport(cl,  varlist=c("time",all.vars(model$call),
                   all.vars(model$drift)[which(all.vars(model$drift)!="x" & all.vars(model$drift)!="t")],
                   all.vars(model$diffusion)[which(all.vars(model$diffusion)!="x" & all.vars(model$diffusion)!="t")]),envir = environment())
                if (base::RNGkind()[1L] == "L'Ecuyer-CMRG") parallel::clusterSetRNGStream(cl)
                rand <- do.call("rbind", parallel::parLapply(cl, 1:R, run1d))
                parallel::stopCluster(cl)
            }else{
                parallel::clusterEvalQ(cl, library(Sim.DiffProc))
                parallel::clusterExport(cl, varlist=c("time",all.vars(model$call),
                   all.vars(model$drift)[which(all.vars(model$drift)!="x" & all.vars(model$drift)!="t")],
                   all.vars(model$diffusion)[which(all.vars(model$diffusion)!="x" & all.vars(model$diffusion)!="t")]),envir = environment())
                if (base::RNGkind()[1L] == "L'Ecuyer-CMRG") parallel::clusterSetRNGStream(cl)
                rand <- do.call("rbind", parallel::parLapply(cl, 1:R, run1d))
                parallel::stopCluster(cl)
             }
         }}else{
         for (i in 1:R){
         rand[i,] = run1d()
         }
        }
        Stat      <- do.call("cbind",lapply(1:R,function(i) fn1d(rand,i,...)))
      }else if (class(model)=="snssde2d"){
       statistic <- match.fun(statistic)
       fn2d      <- function(x,...) statistic(x, ...)
       rand      <- matrix(0,nrow=2*R,ncol=model$M)
       run2d     <- function(...) t(Sim.DiffProc::rsde2d(eval(model$call),at=time))
       if (ncpus > 1L && (have_mc || have_snow)){
         if (have_mc) {
		   base::RNGkind("L'Ecuyer-CMRG")
           rand <- do.call("rbind",parallel::mclapply(1:R,run2d,mc.cores = ncpus,mc.set.seed = TRUE))
         }else if (have_snow) {
            list(...)
            if (is.null(cl)==TRUE) {
                cl  <- parallel::makePSOCKcluster(rep("localhost", ncpus))
                parallel::clusterEvalQ(cl, library(Sim.DiffProc))
                parallel::clusterExport(cl,  varlist=c("time",all.vars(model$call),
                    all.vars(model$driftx)[which(all.vars(model$driftx)!="x" & all.vars(model$driftx)!="y" & all.vars(model$driftx)!="t")],
                    all.vars(model$drifty)[which(all.vars(model$drifty)!="x" & all.vars(model$drifty)!="y" & all.vars(model$drifty)!="t")],
                    all.vars(model$diffx)[which(all.vars(model$diffx)!="x" & all.vars(model$diffx)!="y" & all.vars(model$diffx)!="t")],
                    all.vars(model$diffy)[which(all.vars(model$diffy)!="x" & all.vars(model$diffy)!="y" & all.vars(model$diffy)!="t")]),envir = environment())
                if (base::RNGkind()[1L] == "L'Ecuyer-CMRG") parallel::clusterSetRNGStream(cl)
                rand <- do.call("rbind", parallel::parLapply(cl, 1:R, run2d))
                parallel::stopCluster(cl)
            }else{
                parallel::clusterEvalQ(cl, library(Sim.DiffProc))
                parallel::clusterExport(cl,  varlist=c("time",all.vars(model$call),
                    all.vars(model$driftx)[which(all.vars(model$driftx)!="x" & all.vars(model$driftx)!="y" & all.vars(model$driftx)!="t")],
                    all.vars(model$drifty)[which(all.vars(model$drifty)!="x" & all.vars(model$drifty)!="y" & all.vars(model$drifty)!="t")],
                    all.vars(model$diffx)[which(all.vars(model$diffx)!="x" & all.vars(model$diffx)!="y" & all.vars(model$diffx)!="t")],
                    all.vars(model$diffy)[which(all.vars(model$diffy)!="x" & all.vars(model$diffy)!="y" & all.vars(model$diffy)!="t")]),envir = environment())
                if (base::RNGkind()[1L] == "L'Ecuyer-CMRG") parallel::clusterSetRNGStream(cl)
                rand <- do.call("rbind", parallel::parLapply(cl, 1:R, run2d))
                parallel::stopCluster(cl)
             }
       }}else{
       for (i in 1:R){
          rand[c(seq(2*i-1,2*i,by=1)),]=run2d()
         }
       }
        rand1 <- list(x=rand[seq(1,2*R,by=2),],y=rand[seq(2,2*R,by=2),])
        rand2 <- lapply(1:R,function(i) data.frame(x=rand1$x[i,], y=rand1$y[i,]))
        Stat  <- do.call("cbind",lapply(1:R,function(i) fn2d(rand2[[i]],...)))
      }else if (class(model)=="snssde3d"){
        statistic <- match.fun(statistic)
        fn3d      <- function(x,...) statistic(x, ...)
        rand      <- matrix(0,nrow=3*R,ncol=model$M)
        run3d     <- function(...) t(Sim.DiffProc::rsde3d(eval(model$call),at=time))
        if (ncpus > 1L && (have_mc || have_snow)){
         if (have_mc) {
		   base::RNGkind("L'Ecuyer-CMRG")
           rand <- do.call("rbind",parallel::mclapply(1:R,run3d,mc.cores = ncpus,mc.set.seed = TRUE))
         }else if (have_snow) {
            list(...)
            if (is.null(cl)==TRUE) {
                cl  <- parallel::makePSOCKcluster(rep("localhost", ncpus))
                parallel::clusterEvalQ(cl, library(Sim.DiffProc))
                parallel::clusterExport(cl,  varlist=c("time",all.vars(model$call),
                    all.vars(model$driftx)[which(all.vars(model$driftx)!="x" & all.vars(model$driftx)!="y" & all.vars(model$driftx)!="z" & all.vars(model$driftx)!="t")],
                    all.vars(model$drifty)[which(all.vars(model$drifty)!="x" & all.vars(model$drifty)!="y" & all.vars(model$drifty)!="z" & all.vars(model$drifty)!="t")],
                    all.vars(model$driftz)[which(all.vars(model$driftz)!="x" & all.vars(model$driftz)!="y" & all.vars(model$driftz)!="z" & all.vars(model$driftz)!="t")],
                    all.vars(model$diffx)[which(all.vars(model$diffx)!="x" & all.vars(model$diffx)!="y" & all.vars(model$diffx)!="z" & all.vars(model$diffx)!="t")],
                    all.vars(model$diffy)[which(all.vars(model$diffy)!="x" & all.vars(model$diffy)!="y" & all.vars(model$diffy)!="z" & all.vars(model$diffy)!="t")],
                    all.vars(model$diffz)[which(all.vars(model$diffz)!="x" & all.vars(model$diffz)!="y" & all.vars(model$diffz)!="z" & all.vars(model$diffz)!="t")]),envir = environment())
                if (base::RNGkind()[1L] == "L'Ecuyer-CMRG") parallel::clusterSetRNGStream(cl)
                rand <- do.call("rbind", parallel::parLapply(cl, 1:R, run3d))
                parallel::stopCluster(cl)
            }else{
                parallel::clusterEvalQ(cl, library(Sim.DiffProc))
                parallel::clusterExport(cl,  varlist=c("time",all.vars(model$call),
                    all.vars(model$driftx)[which(all.vars(model$driftx)!="x" & all.vars(model$driftx)!="y" & all.vars(model$driftx)!="z" & all.vars(model$driftx)!="t")],
                    all.vars(model$drifty)[which(all.vars(model$drifty)!="x" & all.vars(model$drifty)!="y" & all.vars(model$drifty)!="z" & all.vars(model$drifty)!="t")],
                    all.vars(model$driftz)[which(all.vars(model$driftz)!="x" & all.vars(model$driftz)!="y" & all.vars(model$driftz)!="z" & all.vars(model$driftz)!="t")],
                    all.vars(model$diffx)[which(all.vars(model$diffx)!="x" & all.vars(model$diffx)!="y" & all.vars(model$diffx)!="z" & all.vars(model$diffx)!="t")],
                    all.vars(model$diffy)[which(all.vars(model$diffy)!="x" & all.vars(model$diffy)!="y" & all.vars(model$diffy)!="z" & all.vars(model$diffy)!="t")],
                    all.vars(model$diffz)[which(all.vars(model$diffz)!="x" & all.vars(model$diffz)!="y" & all.vars(model$diffz)!="z" & all.vars(model$diffz)!="t")]),envir = environment())
                if (base::RNGkind()[1L] == "L'Ecuyer-CMRG") parallel::clusterSetRNGStream(cl)
                rand <- do.call("rbind", parallel::parLapply(cl, 1:R, run3d))
                parallel::stopCluster(cl)
             }
       }}else{
        for (i in 1:R){
          rand[c(seq(3*i-2,3*i,by=1)),]=run3d()
          }
		}
        rand1 <- list(x=rand[seq(1,3*R,by=3),],y=rand[seq(2,3*R,by=3),],z=rand[seq(3,3*R,by=3),])
        rand3 <- lapply(1:R,function(i) data.frame(x=rand1$x[i,], y=rand1$y[i,],z=rand1$z[i,]))
        Stat  <- do.call("cbind",lapply(1:R,function(i) fn3d(rand3[[i]],...)))
       }
       Est    <- round(apply(Stat,1,mean),digits=options()$digits)
       #Sdev  <- round(apply(Stat,1,sd),options()$digits)
       SErr   <- round(apply(Stat,1,sd)/sqrt(R),digits=options()$digits)
       INF    <- round(Est - qnorm(1-0.5*(1-level))*SErr,digits=options()$digits)
       SUP    <- round(Est + qnorm(1-0.5*(1-level))*SErr,digits=options()$digits)
       Conf   <- paste("(",INF,",",SUP,")",sep=" ")
       if (!is.null(exact)){
	      rmse_f <- function(error) sqrt(mean(error^2))
          Exact = round(as.numeric(exact),digits=options()$digits)
          Bias = round(Est - Exact,digits=options()$digits)
          Rmse = round(apply(Stat-Exact,1, rmse_f ),digits=options()$digits)
          TAB  <- data.frame(Exact,Est,Bias,SErr,Rmse,Conf)
          names(TAB) <- c("Exact","Estimate","Bias","Std.Error","RMSE",paste("CI(",100*(1 - level)/2,"%",",",100-100*(1 - level)/2,"%",")",sep=" "))
          if (!is.null(names(exact)) && is.null(names)) {rownames(TAB) <- names(exact)} 
          else if (is.null(names(exact)) && !is.null(names)) {rownames(TAB) <- names} 
          else if (!is.null(names(exact)) && !is.null(names)) {rownames(TAB) <- names} 
          else{rownames(TAB) <- paste("mu",1:length(Est),sep="")}
       }else{
          TAB  <- data.frame(Est,SErr,Conf)
          names(TAB) <- c("Estimate","Std.Error",paste("CI(",100*(1 - level)/2,"%",",",100-100*(1 - level)/2,"%",")",sep=" "))
          if (!is.null(names)) {rownames(TAB) <- names}else{rownames(TAB) <- paste("mu",1:length(Est),sep="")}
       }
structure(list(MC=TAB,name=dimnames(TAB)[[1]],ech=Stat,mod=model,Fn=statistic,corrmat=model$corrmat,time=time,call=match.call(),infC=INF,supC=SUP,dim=model$dim,Class=class(model)),class="MCM.sde")
}


print.MCM.sde <- function(x, digits=NULL, ...)
          {
    class(x) <- "MCM.sde"
	Ito = "It\xf4"
    Encoding(Ito) <- "latin1"
    if (is.null(x$corrmat)){ 
    if (x$Class=="snssde1d"){
    Dr <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)))), list(e = x$mod$drift))))
    DD <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)))), list(e = x$mod$diffusion))))
    if(x$mod$type=="ito"){
    cat(Ito," Sde 1D:","\n",
        " | dX(t) = ", Dr," * dt + ", DD," * dW(t)","\n",
		" | t in [",format(x$mod$t0,digits=digits),",",format(x$mod$T,digits=digits),"] with mesh equal to ",format(x$mod$Dt,digits=digits),"\n",
        sep="")}else{
    cat("Stratonovich Sde 1D:","\n",
        " | dX(t) = ", Dr," * dt + ", DD," o dW(t)","\n",
		" | t in [",format(x$mod$t0,digits=digits),",",format(x$mod$T,digits=digits),"] with mesh equal to ",format(x$mod$Dt,digits=digits),"\n",
        sep="")}
    }else if (x$Class=="snssde2d"){
	Drx <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)))), list(e = x$mod$driftx))))
    DDx <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)))), list(e = x$mod$diffx))))
	Dry <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)))), list(e = x$mod$drifty))))
    DDy <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)))), list(e = x$mod$diffy))))
        if(x$mod$type=="ito"){
    cat(Ito," Sde 2D:","\n",
        " | dX(t) = ", Drx," * dt + ", DDx," * dW1(t)","\n",
        " | dY(t) = ", Dry," * dt + ", DDy," * dW2(t)","\n",
		" | t in [",format(x$mod$t0,digits=digits),",",format(x$mod$T,digits=digits),"] with mesh equal to ",format(x$mod$Dt,digits=digits),"\n",
        sep="")}else{
    cat("Stratonovich Sde 2D:","\n",
        " | dX(t) = ", Drx," * dt + ", DDx," o dW1(t)","\n",
        " | dY(t) = ", Dry," * dt + ", DDy," o dW2(t)","\n",
		" | t in [",format(x$mod$t0,digits=digits),",",format(x$mod$T,digits=digits),"] with mesh equal to ",format(x$mod$Dt,digits=digits),"\n",
        sep="")}
    }else if (x$Class=="snssde3d"){
    Drx <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)),z=quote(Z(t)))), list(e = x$mod$driftx))))
    DDx <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)),z=quote(Z(t)))), list(e = x$mod$diffx))))
	Dry <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)),z=quote(Z(t)))), list(e = x$mod$drifty))))
    DDy <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)),z=quote(Z(t)))), list(e = x$mod$diffy))))
    Drz <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)),z=quote(Z(t)))), list(e = x$mod$driftz))))
    DDz <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)),z=quote(Z(t)))), list(e = x$mod$diffz))))
    if(x$mod$type=="ito"){
    cat(Ito," Sde 3D:","\n",
        " | dX(t) = ", Drx," * dt + ", DDx," * dW1(t)","\n",
        " | dY(t) = ", Dry," * dt + ", DDy," * dW2(t)","\n",
        " | dZ(t) = ", Drz," * dt + ", DDz," * dW3(t)","\n",
		" | t in [",format(x$mod$t0,digits=digits),",",format(x$mod$T,digits=digits),"] with mesh equal to ",format(x$mod$Dt,digits=digits),"\n",
        sep="")}else{
    cat("Stratonovich Sde 3D:","\n",
        " | dX(t) = ", Drx," * dt + ", DDx," o dW1(t)","\n",
        " | dY(t) = ", Dry," * dt + ", DDy," o dW2(t)","\n",
        " | dZ(t) = ", Drz," * dt + ", DDz," o dW3(t)","\n",
		" | t in [",format(x$mod$t0,digits=digits),",",format(x$mod$T,digits=digits),"] with mesh equal to ",format(x$mod$Dt,digits=digits),"\n",
        sep="")}
    }} else {
    if (x$Class=="snssde2d"){
	Drx <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)))), list(e = x$mod$driftx))))
    DDx <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)))), list(e = x$mod$diffx))))
	Dry <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)))), list(e = x$mod$drifty))))
    DDy <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)))), list(e = x$mod$diffy))))
        if(x$mod$type=="ito"){
    cat(Ito," Sde 2D:","\n",
        " | dX(t) = ", Drx," * dt + ", DDx," * dB1(t)","\n",
        " | dY(t) = ", Dry," * dt + ", DDy," * dB2(t)","\n",
		" | t in [",format(x$mod$t0,digits=digits),",",format(x$mod$T,digits=digits),"] with mesh equal to ",format(x$mod$Dt,digits=digits),"\n",
		" | Correlation structure:",sep="")		 
		   prmatrix(x$mod$corrmat,rowlab = rep("            ", 2), collab = rep("", 2),digits=digits)
     }else{
    cat("Stratonovich Sde 2D:","\n",
        " | dX(t) = ", Drx," * dt + ", DDx," o dB1(t)","\n",
        " | dY(t) = ", Dry," * dt + ", DDy," o dB2(t)","\n",
		" | t in [",format(x$mod$t0,digits=digits),",",format(x$mod$T,digits=digits),"] with mesh equal to ",format(x$mod$Dt,digits=digits),"\n",
		" | Correlation structure:",sep="")		 
		   prmatrix(x$mod$corrmat,rowlab = rep("            ", 2), collab = rep("", 2),digits=digits)
		   }
    }else if (x$Class=="snssde3d"){
    Drx <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)),z=quote(Z(t)))), list(e = x$mod$driftx))))
    DDx <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)),z=quote(Z(t)))), list(e = x$mod$diffx))))
	Dry <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)),z=quote(Z(t)))), list(e = x$mod$drifty))))
    DDy <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)),z=quote(Z(t)))), list(e = x$mod$diffy))))
    Drz <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)),z=quote(Z(t)))), list(e = x$mod$driftz))))
    DDz <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)),z=quote(Z(t)))), list(e = x$mod$diffz))))
    if(x$mod$type=="ito"){
    cat(Ito," Sde 3D:","\n",
        " | dX(t) = ", Drx," * dt + ", DDx," * dB1(t)","\n",
        " | dY(t) = ", Dry," * dt + ", DDy," * dB2(t)","\n",
        " | dZ(t) = ", Drz," * dt + ", DDz," * dB3(t)","\n",
		" | t in [",format(x$mod$t0,digits=digits),",",format(x$mod$T,digits=digits),"] with mesh equal to ",format(x$mod$Dt,digits=digits),"\n",
		" | Correlation structure:",sep="")		 
		   prmatrix(x$mod$corrmat,rowlab = rep("      ", 3), collab = rep("", 3),digits=digits)
		}else{
    cat("Stratonovich Sde 3D:","\n",
        " | dX(t) = ", Drx," * dt + ", DDx," o dB1(t)","\n",
        " | dY(t) = ", Dry," * dt + ", DDy," o dB2(t)","\n",
        " | dZ(t) = ", Drz," * dt + ", DDz," o dB3(t)","\n",
		" | t in [",format(x$mod$t0,digits=digits),",",format(x$mod$T,digits=digits),"] with mesh equal to ",format(x$mod$Dt,digits=digits),"\n",
		" | Correlation structure:",sep="")		 
		   prmatrix(x$mod$corrmat,rowlab = rep("      ", 3), collab = rep("", 3),digits=digits)
		   }
    }	
	}
    cat("\nPMCM Based on ", dim(x$ech)[2]," Batches with ",x$mod$M, "-Realisations at time ",x$time,":","\n\n",sep="")
    print(x$MC)
}


.plot.MCM.sde <- function(x,index = 1,type=c("all","hist","qqplot","boxplot","CI"),...)
         {
    class(x) <- "MCM.sde"
    data <- x$ech
    name <- x$name
    if (length(index) > 1) {
       index <- index[1]
       warning("The first element of 'index' will be used")}
    if (index >  dim(data)[1] | index < 1 ) 
	        stop("Subscript out of bounds")
    #if (is.null(x$call$names)) {name <- paste("t",1:dim(data)[1],"*",sep="")}else{
    #name <- sapply(2:length(x$call$names),function(i) x$call$names[[i]])}
    type <- match.arg(type)
    if (type=="all"){par(mfrow=c(2,2))}
    if (any(type == "all" | type == "hist")){
       MASS::truehist(data[index,],xlab = name[index],col="gold",ylab="Density",bty="o",las=1,...)
       abline(v=mean(data[index,]),lty=2,col="red3")
       }
    if (any(type == "all" | type == "qqplot")){
       qq <- qnorm((seq_len(dim(data)[2]))/(dim(data)[2]+1))
       qlab <-"Quantiles of Standard Normal"
       stats::qqplot(qq,data[index,],xlab=qlab,ylab=name[index],las=1)
       abline(a=mean(data[index,]),b=sqrt(var(data[index,])),lty=2,col="red3",lwd=2)
      }
    if (any(type == "all" | type == "boxplot")){
       bx <- boxplot(data[index,],plot=FALSE,names=name[index])
       graphics::bxp(bx, notch = FALSE,las=1, boxfill = "gold", frame.plot = TRUE,outline = FALSE,show.names = TRUE)
      }
    if (any(type == "all" | type == "CI")){
       plot(mean(data[index,]),type="n" ,xaxt='n',yaxt="n",ylim=c(x$infC[index], x$supC[index]),ylab='Confidence Interval',xlab="",las=1,panel.first = grid())
       axis(1, at=1, labels=name[index],lty = 2, lwd = 0.5)
       axis(side=2,lwd=0.5,las=1,cex=0.5)
       arrows(1,x$infC[index],1,x$supC[index],code=3,length=0.25,lwd=2,angle=90,col='red3')
       points(1,mean(data[index,]),pch=21,col="red3",bg="gold",cex=1.5)
      }
    par(mfrow=c(1,1))
    invisible(x)
}

plot.MCM.sde <- function(x,index = 1,type=c("all","hist","qqplot","boxplot","CI"),...) .plot.MCM.sde(x,index,type,...)

