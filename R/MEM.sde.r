## Mon May 27 03:39:04 2019
## Original file Copyright © 2019 A.C. Guidoum, K. Boukhetala
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


MEM.sde <- function(drift, diffusion, ...)  UseMethod("MEM.sde")

MEM.sde.default <- function(drift, diffusion, type=c("ito","str"), solve = FALSE, parms = NULL, init = NULL, time = NULL,...)
                 {
    if (length(drift) != length(diffusion)) 
	                 stop("coefficient of 'drift' and 'diffusion' have a different dimensions")
    if (length(drift) == 1 && length(diffusion) == 1){
         if (!is.expression(drift) || !is.expression(diffusion)) 
		             stop(" coefficient of 'drift' and 'diffusion' must be expressions in 't' and 'x'")}
    else if (length(drift) == 2 && length(diffusion) == 2){
         if (!is.expression(drift) || !is.expression(diffusion)) 
		             stop(" coefficient of 'drift' and 'diffusion' must be expressions in 't', 'x' and 'y'")}
    else if (length(drift) == 3 && length(diffusion) == 3){
         if (!is.expression(drift) || !is.expression(diffusion)) 
		             stop(" coefficient of 'drift' and 'diffusion' must be expressions in 't', 'x', 'y' and 'z'")}
    else {
	    stop("coefficient of 'drift' and 'diffusion' have a dimensions > 3d")}
    if (missing(type)) 
	         type <- "ito"
    if (type=="ito"){nu = 1}else{nu = 2}
    if (solve){
        if (!requireNamespace("deSolve", quietly = TRUE)) {
                cat("The deSolve package is not available, you need to install the package.\\n")
                solve <- FALSE }
      }
    if (solve){
    if (is.null(init))
	    stop("argument 'init' is missing, with no default")
    if (is.null(time))
	    stop("argument 'time' is missing, with no default")
	at = tail(time, n=1)
    }
    if (!is.null(parms)){
       allnames <-names(parms)
         for (i in 1:length(parms)) {
		   currname  <- allnames[i]
		   value     <- toString(parms[[i]])
		   drift     <- parse(text=gsub(pattern = currname, replacement = value, x = drift))
           diffusion <- parse(text=gsub(pattern = currname, replacement = value, x = diffusion))
	   }
    }	
    if (length(drift) == 1 && length(diffusion) == 1){
    ## Means
    Eqm01 <- eval(Simplify(substitute(expression(e1 - 0.5 * (nu-1) * e2 * de2), list(e1 = drift[[1]], e2 = diffusion[[1]], de2 = Deriv(diffusion,"x",cache.exp=FALSE)[[1]],nu=nu))))
    Eqm02 <- eval(Simplify(substitute(expression(0.5 * e1 * S), list(e1=Deriv(Eqm01[[1]],"x",nderiv=2,cache.exp=FALSE),S=quote(S)))))
    Eqm03 <- eval(Simplify(substitute(expression(e1+e2),list(e1 = Eqm01[[1]] , e2 = Eqm02[[1]]))))
    Eqm0F <- eval(Simplify(substitute(substitute(e, list(x=quote(m))), list(e = Eqm03[[1]]))))
    ## Variances
    EqS01 <- eval(Simplify(substitute(expression(e1 - 0.5 * (nu-1) * e2 * de2), list(e1 = drift[[1]], e2 = diffusion[[1]], de2 = Deriv(diffusion,"x",cache.exp=FALSE)[[1]],nu=nu))))
    EqS02 <- eval(Simplify(substitute(expression(2*e1*S),list(e1=Deriv(EqS01[[1]],"x",cache.exp=FALSE),S=quote(S)))))
    EqS03 <- eval(Simplify(substitute(expression(e2*e2),list(e2 = diffusion[[1]]))))
    EqS04 <- eval(Simplify(substitute(expression((e2*dde2+de2*de2)*S),list(e2 = diffusion[[1]],de2 = Deriv(diffusion,"x",cache.exp=FALSE)[[1]], dde2 = Deriv(diffusion,"x",nderiv=2,cache.exp=FALSE)[[1]],S=quote(S)))))
    EqS06 <- eval(Simplify(substitute(expression(e1+e2+e3),list(e1 = EqS02[[1]] , e2 = EqS03[[1]], e3 = EqS04[[1]]))))
    EqS0F <- eval(Simplify(substitute(substitute(e, list(x=quote(m))), list(e = EqS06[[1]]))))
    ## Solve ODE
    if (solve==TRUE){
    mod1d <- function(t,State,Pars){
       with(as.list(c(State, Pars)),{
          dm = eval(Eqm0F)
          dS = eval(EqS0F)
          return(list(c(dm,dS)))
          })
    }
    Sol <- deSolve::ode(func = mod1d, y = init,parms = parms, times = time,...)}else{
    Sol <- NULL}
    ## output
    structure(list(Means=Eqm0F,Var=EqS0F,drift=drift,diffusion=diffusion,type=type,dim=length(drift),sol.ode=Sol,call=match.call()),class="MEM.sde")
    }else if (length(drift) == 2 && length(diffusion) == 2){
    ## Means
    # m1(t)
    Eqm101 <- eval(Simplify(substitute(expression(e1 - 0.5 * (nu-1) * e2 * de2), list(e1 = drift[1][[1]], e2 = diffusion[1][[1]], de2 = Deriv(diffusion[1],"x",cache.exp=FALSE)[[1]],nu=nu))))
    Eqm102 <- eval(Simplify(substitute(expression(0.5 * e1 * S1 + 0.5 * e2 * S2 + e3 * C12), 
               list(e1=Deriv(Eqm101[[1]],"x",nderiv=2,cache.exp=FALSE),
                    e2=Deriv(Eqm101[[1]],"y",nderiv=2,cache.exp=FALSE),
                    e3=Deriv(Deriv(Eqm101[[1]],"x",cache.exp=FALSE),"y",cache.exp=FALSE),
                    S1=quote(S1),S2=quote(S2),C12=quote(C12)))))
    Eqm103 <- eval(Simplify(substitute(expression(e1+e2),list(e1 = Eqm101[[1]] , e2 = Eqm102[[1]]))))
    Eqm10F <- eval(Simplify(substitute(substitute(e, list(x=quote(m1),y=quote(m2))), list(e = Eqm103[[1]]))))
    # m2(t)
    Eqm201 <- eval(Simplify(substitute(expression(e1 - 0.5 * (nu-1) * e2 * de2), list(e1 = drift[2][[1]], e2 = diffusion[2][[1]], de2 = Deriv(diffusion[2],"y",cache.exp=FALSE)[[1]],nu=nu))))
    Eqm202 <- eval(Simplify(substitute(expression(0.5 * e1 * S1 + 0.5 * e2 * S2 + e3 * C12), 
               list(e1=Deriv(Eqm201[[1]],"x",nderiv=2,cache.exp=FALSE),
                    e2=Deriv(Eqm201[[1]],"y",nderiv=2,cache.exp=FALSE),
                    e3=Deriv(Deriv(Eqm201[[1]],"x",cache.exp=FALSE),"y",cache.exp=FALSE),
                    S1=quote(S1),S2=quote(S2),C12=quote(C12)))))
    Eqm203 <- eval(Simplify(substitute(expression(e1+e2),list(e1 = Eqm201[[1]] , e2 = Eqm202[[1]]))))
    Eqm20F <- eval(Simplify(substitute(substitute(e, list(x=quote(m1),y=quote(m2))), list(e = Eqm203[[1]]))))
    ## Variances
    # S1(t)
    EqS101 <- eval(Simplify(substitute(expression(e1 - 0.5 * (nu-1) * e2 * de2), list(e1 = drift[1][[1]], e2 = diffusion[1][[1]], de2 = Deriv(diffusion[1],"x",cache.exp=FALSE)[[1]],nu=nu))))
    EqS102 <- eval(Simplify(substitute(expression( ( 2*e1 + e2 * dde2 + de2^2 ) * S1), 
               list(e1   = Deriv(EqS101[[1]],"x",cache.exp=FALSE), 
                    e2   = diffusion[1][[1]], 
                    dde2 = Deriv(diffusion[1],"x",nderiv=2,cache.exp=FALSE)[[1]],
                    de2  = Deriv(diffusion[1],"x",cache.exp=FALSE)[[1]],S1=quote(S1)))))
    EqS103 <- eval(Simplify(substitute(expression( e1 ^2 + ( e2 * dde2 + de2^2 ) * S2), 
               list(e1   = diffusion[1][[1]],
                    e2   = diffusion[1][[1]], 
                    dde2 = Deriv(diffusion[1],"y",nderiv=2,cache.exp=FALSE)[[1]],
                    de2  = Deriv(diffusion[1],"y",cache.exp=FALSE)[[1]],S2=quote(S2)))))
    EqS104 <- eval(Simplify(substitute(expression( 2 * ( e1 + e2 * dde2 + de1 * de2 ) * C12), 
               list(e1   = Deriv(EqS101[[1]],"y",cache.exp=FALSE), 
                    e2   = diffusion[1][[1]], 
                    dde2 = Deriv(Deriv(diffusion[1],"x",cache.exp=FALSE),"y",cache.exp=FALSE)[[1]],
                    de1  = Deriv(diffusion[1],"x",cache.exp=FALSE)[[1]],
                    de2  = Deriv(diffusion[1],"y",cache.exp=FALSE)[[1]],C12=quote(C12)))))
    EqS105 <- eval(Simplify(substitute(expression(e1+e2+e3),list(e1 = EqS102[[1]] , e2 = EqS103[[1]], e3 = EqS104[[1]]))))
    EqS10F <- eval(Simplify(substitute(substitute(e, list(x=quote(m1),y=quote(m2))), list(e = EqS105[[1]]))))
    # S2(t)
    EqS201 <- eval(Simplify(substitute(expression(e1 - 0.5 * (nu-1) * e2 * de2), list(e1 = drift[2][[1]], e2 = diffusion[2][[1]], de2 = Deriv(diffusion[2],"y",cache.exp=FALSE)[[1]],nu=nu))))
    EqS202 <- eval(Simplify(substitute(expression( ( 2*e1 + e2 * dde2 + de2^2 ) * S2), 
               list(e1   = Deriv(EqS201[[1]],"y",cache.exp=FALSE), 
                    e2   = diffusion[2][[1]], 
                    dde2 = Deriv(diffusion[2],"y",nderiv=2,cache.exp=FALSE)[[1]],
                    de2  = Deriv(diffusion[2],"y",cache.exp=FALSE)[[1]],S2=quote(S2)))))
    EqS203 <- eval(Simplify(substitute(expression( e1 ^2 + ( e2 * dde2 + de2^2 ) * S1), 
               list(e1   = diffusion[2][[1]],
                    e2   = diffusion[2][[1]], 
                    dde2 = Deriv(diffusion[2],"x",nderiv=2,cache.exp=FALSE)[[1]],
                    de2  = Deriv(diffusion[2],"x",cache.exp=FALSE)[[1]],S1=quote(S1)))))
    EqS204 <- eval(Simplify(substitute(expression( 2 * ( e1 + e2 * dde2 + de1 * de2 ) * C12), 
               list(e1   = Deriv(EqS201[[1]],"x",cache.exp=FALSE), 
                    e2   = diffusion[2][[1]], 
                    dde2 = Deriv(Deriv(diffusion[2],"y",cache.exp=FALSE),"x",cache.exp=FALSE)[[1]],
                    de1  = Deriv(diffusion[2],"y",cache.exp=FALSE)[[1]],
                    de2  = Deriv(diffusion[2],"x",cache.exp=FALSE)[[1]],C12=quote(C12)))))
    EqS205 <- eval(Simplify(substitute(expression(e1+e2+e3),list(e1 = EqS202[[1]] , e2 = EqS203[[1]], e3 = EqS204[[1]]))))
    EqS20F <- eval(Simplify(substitute(substitute(e, list(x=quote(m1),y=quote(m2))), list(e = EqS205[[1]]))))
    ## Covariance
    # C12(t)
    Z1 <- eval(Simplify(substitute(expression(e1 - 0.5 * (nu-1) * e2 * de2), list(e1 = drift[1][[1]], e2 = diffusion[1][[1]], de2 = Deriv(diffusion[1],"x",cache.exp=FALSE)[[1]],nu=nu))))
    Z2 <- eval(Simplify(substitute(expression(e1 - 0.5 * (nu-1) * e2 * de2), list(e1 = drift[2][[1]], e2 = diffusion[2][[1]], de2 = Deriv(diffusion[2],"y",cache.exp=FALSE)[[1]],nu=nu))))
    EqC1201 <- eval(Simplify(substitute(expression( ( e + 0.5 * e2 * dde1 + de1 * de2  + 0.5 * e1 * dde2) * S1), 
               list(e    = Deriv(Z2[[1]],"x",cache.exp=FALSE), 
                    e2   = diffusion[2][[1]], 
                    dde1 = Deriv(diffusion[1],"x",nderiv=2,cache.exp=FALSE)[[1]],
                    de1  = Deriv(diffusion[1],"x",cache.exp=FALSE)[[1]],
                    de2  = Deriv(diffusion[2],"x",cache.exp=FALSE)[[1]],
                    e1   = diffusion[1][[1]],
                    dde2 = Deriv(diffusion[2],"x",nderiv=2,cache.exp=FALSE)[[1]],
                    S1=quote(S1)))))
    EqC1202 <- eval(Simplify(substitute(expression( ( e + 0.5 * e2 * dde1 + de1 * de2  + 0.5 * e1 * dde2) * S2), 
               list(e    = Deriv(Z1[[1]],"y",cache.exp=FALSE), 
                    e2   = diffusion[2][[1]], 
                    dde1 = Deriv(diffusion[1],"y",nderiv=2,cache.exp=FALSE)[[1]],
                    de1  = Deriv(diffusion[1],"y",cache.exp=FALSE)[[1]],
                    de2  = Deriv(diffusion[2],"y",cache.exp=FALSE)[[1]],
                    e1   = diffusion[1][[1]],
                    dde2 = Deriv(diffusion[2],"y",nderiv=2,cache.exp=FALSE)[[1]],
                    S2=quote(S2)))))
    EqC1203 <- eval(Simplify(substitute(expression( ( e1 + e2 + e4 * dde3 + de3 * de4 + de5 * de6 + e3 * dde4) * C12), 
               list(e1   = Deriv(Z1[[1]],"x",cache.exp=FALSE), 
                    e2   = Deriv(Z2[[1]],"y",cache.exp=FALSE),
                    e3   = diffusion[1][[1]],
                    e4   = diffusion[2][[1]], 
                    dde3 = Deriv(Deriv(diffusion[1],"y",cache.exp=FALSE),"x",cache.exp=FALSE)[[1]],
                    dde4 = Deriv(Deriv(diffusion[2],"y",cache.exp=FALSE),"x",cache.exp=FALSE)[[1]],
                    de3  = Deriv(diffusion[1],"y",cache.exp=FALSE)[[1]],
                    de4  = Deriv(diffusion[2],"x",cache.exp=FALSE)[[1]],
                    de5  = Deriv(diffusion[1],"x",cache.exp=FALSE)[[1]],
                    de6  = Deriv(diffusion[2],"y",cache.exp=FALSE)[[1]],
                    C12=quote(C12)))))
    EqC1204 <- eval(Simplify(substitute(expression(e1 * e2), list(e1 = diffusion[1][[1]], e2 = diffusion[2][[1]]))))
    EqC1205 <- eval(Simplify(substitute(expression(e1 + e2 + e3 + e4), list(e1 = EqC1201[[1]], e2 = EqC1202[[1]],e3 = EqC1203[[1]], e4 = EqC1204[[1]]))))
    EqC120F <- eval(Simplify(substitute(substitute(e, list(x=quote(m1),y=quote(m2))), list(e = EqC1205[[1]]))))
    ## Solve ODE
    if (solve==TRUE){
    mod2d <- function(t,State,Pars){
       with(as.list(c(State, Pars)),{
          dm1  = eval(Eqm10F)
          dm2  = eval(Eqm20F)
          dS1  = eval(EqS10F)
          dS2  = eval(EqS20F)
          dC12 = eval(EqC120F)
          return(list(c(dm1,dm2,dS1,dS2,dC12)))
          })
    }
    Sol <- deSolve::ode(func = mod2d, y = init,parms = parms, times = time,...)
	M = c(paste("| m1(",at,")  = ",sep=""),paste("| m2(",at,")  = ",sep=""))
    S = c(paste("| S1(",at,")  = ",sep=""),paste("| S2(",at,")  = ",sep=""))
    C = c(paste("| C12(",at,")  = ",sep=""),"")
	res1 = round(c(tail(Sol[,"m1"], n=1),tail(Sol[,"m2"], n=1)),options()$digits)
    res2 = round(c(tail(Sol[,"S1"], n=1),tail(Sol[,"S2"], n=1)),options()$digits)
    res3 = c(round(tail(Sol[,"C12"], n=1),options()$digits),"")
	MEMF <- data.frame(M,res1,S,res2,C,res3,fix.empty.names = FALSE,row.names = c(" ", ""))
    }else{
    Sol <- NULL
	MEMF<- NULL}
    ## output
    structure(list(Means=list(m1=Eqm10F,m2=Eqm20F),Var=list(S1=EqS10F,S2=EqS20F,C12=EqC120F),drift=drift,diffusion=diffusion,type=type,dim=length(drift),sol.ode=Sol,TAB=MEMF,call=match.call()),class="MEM.sde")
    }else if (length(drift) == 3 && length(diffusion) == 3){
    ## Means
    # m1(t)
    Eqm101 <- eval(Simplify(substitute(expression(e1 - 0.5 * (nu-1) * e2 * de2), list(e1 = drift[1][[1]], e2 = diffusion[1][[1]], de2 = Deriv(diffusion[1],"x",cache.exp=FALSE)[[1]],nu=nu))))
    Eqm102 <- eval(Simplify(substitute(expression(e1 + 0.5 * (ddxe1 * S1 + ddye1 * S2 + ddze1 * S3) ), 
               list(e1   = Eqm101[[1]],
                    ddxe1= Deriv(Eqm101[[1]],"x",nderiv=2,cache.exp=FALSE),
                    ddye1= Deriv(Eqm101[[1]],"y",nderiv=2,cache.exp=FALSE),
                    ddze1= Deriv(Eqm101[[1]],"z",nderiv=2,cache.exp=FALSE),
                    S1=quote(S1),S2=quote(S2),S3=quote(S3)))))
    Eqm103 <- eval(Simplify(substitute(expression( dxye1 * C12 + dxze1 * C13 + dyze1 * C23 ), 
               list(dxye1= Deriv(Deriv(Eqm101[[1]],"x",cache.exp=FALSE),"y",cache.exp=FALSE),
                    dxze1= Deriv(Deriv(Eqm101[[1]],"x",cache.exp=FALSE),"z",cache.exp=FALSE),
                    dyze1= Deriv(Deriv(Eqm101[[1]],"y",cache.exp=FALSE),"z",cache.exp=FALSE),
                    C12=quote(C12),C13=quote(C13),C23=quote(C23)))))
    Eqm104 <- eval(Simplify(substitute(expression(e1+e2),list(e1 = Eqm102[[1]] , e2 = Eqm103[[1]]))))
    Eqm10F <- eval(Simplify(substitute(substitute(e, list(x=quote(m1),y=quote(m2),z=quote(m3))), list(e = Eqm104[[1]]))))
    # m2(t)
    Eqm201 <- eval(Simplify(substitute(expression(e1 - 0.5 * (nu-1) * e2 * de2), list(e1 = drift[2][[1]], e2 = diffusion[2][[1]], de2 = Deriv(diffusion[2],"y",cache.exp=FALSE)[[1]],nu=nu))))
    Eqm202 <- eval(Simplify(substitute(expression(e1 + 0.5 * (ddxe1 * S1 + ddye1 * S2 + ddze1 * S3) ), 
               list(e1   = Eqm201[[1]],
                    ddxe1= Deriv(Eqm201[[1]],"x",nderiv=2,cache.exp=FALSE),
                    ddye1= Deriv(Eqm201[[1]],"y",nderiv=2,cache.exp=FALSE),
                    ddze1= Deriv(Eqm201[[1]],"z",nderiv=2,cache.exp=FALSE),
                    S1=quote(S1),S2=quote(S2),S3=quote(S3)))))
    Eqm203 <- eval(Simplify(substitute(expression( dxye1 * C12 + dxze1 * C13 + dyze1 * C23 ), 
               list(dxye1= Deriv(Deriv(Eqm201[[1]],"x",cache.exp=FALSE),"y",cache.exp=FALSE),
                    dxze1= Deriv(Deriv(Eqm201[[1]],"x",cache.exp=FALSE),"z",cache.exp=FALSE),
                    dyze1= Deriv(Deriv(Eqm201[[1]],"y",cache.exp=FALSE),"z",cache.exp=FALSE),
                    C12=quote(C12),C13=quote(C13),C23=quote(C23)))))
    Eqm204 <- eval(Simplify(substitute(expression(e1+e2),list(e1 = Eqm202[[1]] , e2 = Eqm203[[1]]))))
    Eqm20F <- eval(Simplify(substitute(substitute(e, list(x=quote(m1),y=quote(m2),z=quote(m3))), list(e = Eqm204[[1]]))))
    # m3(t)
    Eqm301 <- eval(Simplify(substitute(expression(e1 - 0.5 * (nu-1) * e2 * de2), list(e1 = drift[3][[1]], e2 = diffusion[3][[1]], de2 = Deriv(diffusion[3],"z",cache.exp=FALSE)[[1]],nu=nu))))
    Eqm302 <- eval(Simplify(substitute(expression(e1 + 0.5 * (ddxe1 * S1 + ddye1 * S2 + ddze1 * S3) ), 
               list(e1   = Eqm301[[1]],
                    ddxe1= Deriv(Eqm301[[1]],"x",nderiv=2,cache.exp=FALSE),
                    ddye1= Deriv(Eqm301[[1]],"y",nderiv=2,cache.exp=FALSE),
                    ddze1= Deriv(Eqm301[[1]],"z",nderiv=2,cache.exp=FALSE),
                    S1=quote(S1),S2=quote(S2),S3=quote(S3)))))
    Eqm303 <- eval(Simplify(substitute(expression( dxye1 * C12 + dxze1 * C13 + dyze1 * C23 ), 
               list(dxye1= Deriv(Deriv(Eqm301[[1]],"x",cache.exp=FALSE),"y",cache.exp=FALSE),
                    dxze1= Deriv(Deriv(Eqm301[[1]],"x",cache.exp=FALSE),"z",cache.exp=FALSE),
                    dyze1= Deriv(Deriv(Eqm301[[1]],"y",cache.exp=FALSE),"z",cache.exp=FALSE),
                    C12=quote(C12),C13=quote(C13),C23=quote(C23)))))
    Eqm304 <- eval(Simplify(substitute(expression(e1+e2),list(e1 = Eqm302[[1]] , e2 = Eqm303[[1]]))))
    Eqm30F <- eval(Simplify(substitute(substitute(e, list(x=quote(m1),y=quote(m2),z=quote(m3))), list(e = Eqm304[[1]]))))
    ## Variances
    # S1(t)
    EqS101 <- eval(Simplify(substitute(expression(e1 - 0.5 * (nu-1) * e2 * de2), list(e1 = drift[1][[1]], e2 = diffusion[1][[1]], de2 = Deriv(diffusion[1],"x",cache.exp=FALSE)[[1]],nu=nu))))
    EqS102 <- eval(Simplify(substitute(expression( 2 * dxe1 * S1 + 2* dye1 * C12 + 2* dze1 * C13 + e2^2 ), 
               list(dxe1 = Deriv(EqS101[[1]],"x",cache.exp=FALSE), 
                    dye1 = Deriv(EqS101[[1]],"y",cache.exp=FALSE), 
                    dze1 = Deriv(EqS101[[1]],"z",cache.exp=FALSE),
                    e2   = diffusion[1][[1]],
                    S1=quote(S1),C12=quote(C12),C13=quote(C13)))))
    EqS103 <- eval(Simplify(substitute(expression( (e1*dxxe1 + dxe1^2)*S1 + (e1*dyye1+dye1^2)*S2 +(e1*dzze1+dze1^2)*S3), 
               list(e1    = diffusion[1][[1]],
                    dxe1  = Deriv(diffusion[1][[1]],"x",cache.exp=FALSE), 
                    dye1  = Deriv(diffusion[1][[1]],"y",cache.exp=FALSE),
                    dze1  = Deriv(diffusion[1][[1]],"z",cache.exp=FALSE),
                    dxxe1 = Deriv(diffusion[1][[1]],"x",nderiv=2,cache.exp=FALSE),
                    dyye1 = Deriv(diffusion[1][[1]],"y",nderiv=2,cache.exp=FALSE),
                    dzze1 = Deriv(diffusion[1][[1]],"z",nderiv=2,cache.exp=FALSE),
                    S1=quote(S1),S2=quote(S2),S3=quote(S3)))))
    EqS104 <- eval(Simplify(substitute(expression( (2*e1*dxye1+2*dxe1*dye1)*C12 + (2*e1*dxze1+2*dxe1*dze1)*C13 + (2*e1*dyze1+2*dye1*dze1)*C23 ), 
               list(e1    = diffusion[1][[1]],
                    dxe1  = Deriv(diffusion[1][[1]],"x",cache.exp=FALSE), 
                    dye1  = Deriv(diffusion[1][[1]],"y",cache.exp=FALSE),
                    dze1  = Deriv(diffusion[1][[1]],"z",cache.exp=FALSE),
                    dxye1 = Deriv(Deriv(diffusion[1][[1]],"x",cache.exp=FALSE),"y",cache.exp=FALSE),
                    dxze1 = Deriv(Deriv(diffusion[1][[1]],"x",cache.exp=FALSE),"z",cache.exp=FALSE),
                    dyze1 = Deriv(Deriv(diffusion[1][[1]],"y",cache.exp=FALSE),"z",cache.exp=FALSE),
                    C12=quote(C12),C13=quote(C13),C23=quote(C23)))))
    EqS105 <- eval(Simplify(substitute(expression(e1+e2+e3),list(e1 = EqS102[[1]] , e2 = EqS103[[1]], e3 = EqS104[[1]]))))
    EqS10F <- eval(Simplify(substitute(substitute(e, list(x=quote(m1),y=quote(m2),z=quote(m3))), list(e = EqS105[[1]]))))
    # S2(t)
    EqS201 <- eval(Simplify(substitute(expression(e1 - 0.5 * (nu-1) * e2 * de2), list(e1 = drift[2][[1]], e2 = diffusion[2][[1]], de2 = Deriv(diffusion[2],"y",cache.exp=FALSE)[[1]],nu=nu))))
    EqS202 <- eval(Simplify(substitute(expression( 2 * dxe1 * C12 + 2* dye1 * S2 + 2* dze1 * C23 + e2^2 ), 
               list(dxe1 = Deriv(EqS201[[1]],"x",cache.exp=FALSE), 
                    dye1 = Deriv(EqS201[[1]],"y",cache.exp=FALSE), 
                    dze1 = Deriv(EqS201[[1]],"z",cache.exp=FALSE),
                    e2   = diffusion[2][[1]],
                    S2=quote(S2),C12=quote(C12),C23=quote(C23)))))
    EqS203 <- eval(Simplify(substitute(expression( (e1*dxxe1 + dxe1^2)*S1 + (e1*dyye1+dye1^2)*S2 +(e1*dzze1+dze1^2)*S3), 
               list(e1    = diffusion[2][[1]],
                    dxe1  = Deriv(diffusion[2][[1]],"x",cache.exp=FALSE), 
                    dye1  = Deriv(diffusion[2][[1]],"y",cache.exp=FALSE),
                    dze1  = Deriv(diffusion[2][[1]],"z",cache.exp=FALSE),
                    dxxe1 = Deriv(diffusion[2][[1]],"x",nderiv=2,cache.exp=FALSE),
                    dyye1 = Deriv(diffusion[2][[1]],"y",nderiv=2,cache.exp=FALSE),
                    dzze1 = Deriv(diffusion[2][[1]],"z",nderiv=2,cache.exp=FALSE),
                    S1=quote(S1),S2=quote(S2),S3=quote(S3)))))
    EqS204 <- eval(Simplify(substitute(expression( (2*e1*dxye1+2*dxe1*dye1)*C12 + (2*e1*dxze1+2*dxe1*dze1)*C13 + (2*e1*dyze1+2*dye1*dze1)*C23 ), 
               list(e1    = diffusion[2][[1]],
                    dxe1  = Deriv(diffusion[2][[1]],"x",cache.exp=FALSE), 
                    dye1  = Deriv(diffusion[2][[1]],"y",cache.exp=FALSE),
                    dze1  = Deriv(diffusion[2][[1]],"z",cache.exp=FALSE),
                    dxye1 = Deriv(Deriv(diffusion[2][[1]],"x",cache.exp=FALSE),"y",cache.exp=FALSE),
                    dxze1 = Deriv(Deriv(diffusion[2][[1]],"x",cache.exp=FALSE),"z",cache.exp=FALSE),
                    dyze1 = Deriv(Deriv(diffusion[2][[1]],"y",cache.exp=FALSE),"z",cache.exp=FALSE),
                    C12=quote(C12),C13=quote(C13),C23=quote(C23)))))
    EqS205 <- eval(Simplify(substitute(expression(e1+e2+e3),list(e1 = EqS202[[1]] , e2 = EqS203[[1]], e3 = EqS204[[1]]))))
    EqS20F <- eval(Simplify(substitute(substitute(e, list(x=quote(m1),y=quote(m2),z=quote(m3))), list(e = EqS205[[1]]))))
    # S3(t)
    EqS301 <- eval(Simplify(substitute(expression(e1 - 0.5 * (nu-1) * e2 * de2), list(e1 = drift[3][[1]], e2 = diffusion[3][[1]], de2 = Deriv(diffusion[3],"z",cache.exp=FALSE)[[1]],nu=nu))))
    EqS302 <- eval(Simplify(substitute(expression( 2 * dxe1 * C13 + 2* dye1 * C23 + 2* dze1 * S3 + e2^2 ), 
               list(dxe1 = Deriv(EqS301[[1]],"x",cache.exp=FALSE), 
                    dye1 = Deriv(EqS301[[1]],"y",cache.exp=FALSE), 
                    dze1 = Deriv(EqS301[[1]],"z",cache.exp=FALSE),
                    e2   = diffusion[3][[1]],
                    S3=quote(S3),C13=quote(C13),C23=quote(C23)))))
    EqS303 <- eval(Simplify(substitute(expression( (e1*dxxe1 + dxe1^2)*S1 + (e1*dyye1+dye1^2)*S2 +(e1*dzze1+dze1^2)*S3), 
               list(e1    = diffusion[3][[1]],
                    dxe1  = Deriv(diffusion[3][[1]],"x",cache.exp=FALSE), 
                    dye1  = Deriv(diffusion[3][[1]],"y",cache.exp=FALSE),
                    dze1  = Deriv(diffusion[3][[1]],"z",cache.exp=FALSE),
                    dxxe1 = Deriv(diffusion[3][[1]],"x",nderiv=2,cache.exp=FALSE),
                    dyye1 = Deriv(diffusion[3][[1]],"y",nderiv=2,cache.exp=FALSE),
                    dzze1 = Deriv(diffusion[3][[1]],"z",nderiv=2,cache.exp=FALSE),
                    S1=quote(S1),S2=quote(S2),S3=quote(S3)))))
    EqS304 <- eval(Simplify(substitute(expression( (2*e1*dxye1+2*dxe1*dye1)*C12 + (2*e1*dxze1+2*dxe1*dze1)*C13 + (2*e1*dyze1+2*dye1*dze1)*C23 ), 
               list(e1    = diffusion[3][[1]],
                    dxe1  = Deriv(diffusion[3][[1]],"x",cache.exp=FALSE), 
                    dye1  = Deriv(diffusion[3][[1]],"y",cache.exp=FALSE),
                    dze1  = Deriv(diffusion[3][[1]],"z",cache.exp=FALSE),
                    dxye1 = Deriv(Deriv(diffusion[3][[1]],"x",cache.exp=FALSE),"y",cache.exp=FALSE),
                    dxze1 = Deriv(Deriv(diffusion[3][[1]],"x",cache.exp=FALSE),"z",cache.exp=FALSE),
                    dyze1 = Deriv(Deriv(diffusion[3][[1]],"y",cache.exp=FALSE),"z",cache.exp=FALSE),
                    C12=quote(C12),C13=quote(C13),C23=quote(C23)))))
    EqS305 <- eval(Simplify(substitute(expression(e1+e2+e3),list(e1 = EqS302[[1]] , e2 = EqS303[[1]], e3 = EqS304[[1]]))))
    EqS30F <- eval(Simplify(substitute(substitute(e, list(x=quote(m1),y=quote(m2),z=quote(m3))), list(e = EqS305[[1]]))))
    ## Covariance
    # C12(t)
    Z1 <- eval(Simplify(substitute(expression(e1 - 0.5 * (nu-1) * e2 * de2), list(e1 = drift[1][[1]], e2 = diffusion[1][[1]], de2 = Deriv(diffusion[1],"x",cache.exp=FALSE)[[1]],nu=nu))))
    Z2 <- eval(Simplify(substitute(expression(e1 - 0.5 * (nu-1) * e2 * de2), list(e1 = drift[2][[1]], e2 = diffusion[2][[1]], de2 = Deriv(diffusion[2],"y",cache.exp=FALSE)[[1]],nu=nu))))
    EqC1201 <- eval(Simplify(substitute(expression( dxe2 * S1 + dye1 * S2 + ( dxe1+dye2 ) *C12 + dze1 * C23 + dze2 * C13 + g1*g2), 
               list(dxe1 = Deriv(Z1[[1]],"x",cache.exp=FALSE),
                    dye1 = Deriv(Z1[[1]],"y",cache.exp=FALSE), 
                    dze1 = Deriv(Z1[[1]],"z",cache.exp=FALSE), 
                    dxe2 = Deriv(Z2[[1]],"x",cache.exp=FALSE), 
                    dye2 = Deriv(Z2[[1]],"y",cache.exp=FALSE), 
                    dze2 = Deriv(Z2[[1]],"z",cache.exp=FALSE), 
                    g1   = diffusion[1][[1]],
                    g2   = diffusion[2][[1]],  
                    C12=quote(C12),C13=quote(C13),C23=quote(C23),S1=quote(S1),S2=quote(S2)))))
    EqC1202 <- eval(Simplify(substitute(expression( 0.5 * ( g2 * dxxg1 + 2 * dxg1 * dxg2 + g1 * dxxg2)*S1 + 0.5 * ( g2 * dyyg1 + 2 * dyg1 * dyg2 + g1 * dyyg2)*S2 + 0.5 * ( g2 * dzzg1 + 2 * dzg1 * dzg2 + g1 * dzzg2)*S3 ), 
               list(g1   = diffusion[1][[1]],
                    g2   = diffusion[2][[1]], 
                    dxg1 = Deriv(diffusion[1],"x",cache.exp=FALSE)[[1]],
                    dxxg1= Deriv(diffusion[1],"x",nderiv=2,cache.exp=FALSE)[[1]],
                    dyg1 = Deriv(diffusion[1],"y",cache.exp=FALSE)[[1]],
                    dyyg1= Deriv(diffusion[1],"y",nderiv=2,cache.exp=FALSE)[[1]],
                    dzg1 = Deriv(diffusion[1],"z",cache.exp=FALSE)[[1]],
                    dzzg1= Deriv(diffusion[1],"z",nderiv=2,cache.exp=FALSE)[[1]],
                    dxg2 = Deriv(diffusion[2],"x",cache.exp=FALSE)[[1]],
                    dxxg2= Deriv(diffusion[2],"x",nderiv=2,cache.exp=FALSE)[[1]],
                    dyg2 = Deriv(diffusion[2],"y",cache.exp=FALSE)[[1]],
                    dyyg2= Deriv(diffusion[2],"y",nderiv=2,cache.exp=FALSE)[[1]],
                    dzg2 = Deriv(diffusion[2],"z",cache.exp=FALSE)[[1]],
                    dzzg2= Deriv(diffusion[2],"z",nderiv=2,cache.exp=FALSE)[[1]],
                    S1=quote(S1),S2=quote(S2),S3=quote(S3)))))
    EqC1203 <- eval(Simplify(substitute(expression( (g2*dxyg1+dyg1*dxg2+dxg1*dyg2+g1*dxyg2)*C12 + (g2*dxzg1+dzg1*dxg2+dxg1*dzg2+g1*dxzg2)*C13 + (g2*dyzg1+dzg1*dyg2+dyg1*dzg2+g1*dyzg2)*C23 ), 
               list(g1   = diffusion[1][[1]],
                    g2   = diffusion[2][[1]],
                    dxg1 = Deriv(diffusion[1],"x",cache.exp=FALSE)[[1]],
                    dyg1 = Deriv(diffusion[1],"y",cache.exp=FALSE)[[1]],
                    dzg1 = Deriv(diffusion[1],"z",cache.exp=FALSE)[[1]], 
                    dxg2 = Deriv(diffusion[2],"x",cache.exp=FALSE)[[1]],
                    dyg2 = Deriv(diffusion[2],"y",cache.exp=FALSE)[[1]],
                    dzg2 = Deriv(diffusion[2],"z",cache.exp=FALSE)[[1]],
                    dxyg1= Deriv(Deriv(diffusion[1],"x",cache.exp=FALSE),"y",cache.exp=FALSE)[[1]],
                    dxzg1= Deriv(Deriv(diffusion[1],"x",cache.exp=FALSE),"z",cache.exp=FALSE)[[1]],
                    dyzg1= Deriv(Deriv(diffusion[1],"y",cache.exp=FALSE),"z",cache.exp=FALSE)[[1]],
                    dxyg2= Deriv(Deriv(diffusion[2],"x",cache.exp=FALSE),"y",cache.exp=FALSE)[[1]],
                    dxzg2= Deriv(Deriv(diffusion[2],"x",cache.exp=FALSE),"z",cache.exp=FALSE)[[1]],
                    dyzg2= Deriv(Deriv(diffusion[2],"y",cache.exp=FALSE),"z",cache.exp=FALSE)[[1]],
                    C12=quote(C12),C13=quote(C13),C23=quote(C23)))))
    EqC1204 <- eval(Simplify(substitute(expression(e1 + e2 + e3 ), list(e1 = EqC1201[[1]], e2 = EqC1202[[1]],e3 = EqC1203[[1]]))))
    EqC120F <- eval(Simplify(substitute(substitute(e, list(x=quote(m1),y=quote(m2),z=quote(m3))), list(e = EqC1204[[1]]))))
    # C13(t)
    Z3 <- eval(Simplify(substitute(expression(e1 - 0.5 * (nu-1) * e2 * de2), list(e1 = drift[1][[1]], e2 = diffusion[1][[1]], de2 = Deriv(diffusion[1],"x",cache.exp=FALSE)[[1]],nu=nu))))
    Z4 <- eval(Simplify(substitute(expression(e1 - 0.5 * (nu-1) * e2 * de2), list(e1 = drift[3][[1]], e2 = diffusion[3][[1]], de2 = Deriv(diffusion[3],"z",cache.exp=FALSE)[[1]],nu=nu))))
    EqC1301 <- eval(Simplify(substitute(expression( dxe3 * S1 + dze1 * S3 + ( dxe1+dze3 ) *C13 + dye1 * C23 + dye3 * C12 + g1*g3), 
               list(dxe1 = Deriv(Z3[[1]],"x",cache.exp=FALSE),
                    dye1 = Deriv(Z3[[1]],"y",cache.exp=FALSE), 
                    dze1 = Deriv(Z3[[1]],"z",cache.exp=FALSE), 
                    dxe3 = Deriv(Z4[[1]],"x",cache.exp=FALSE), 
                    dye3 = Deriv(Z4[[1]],"y",cache.exp=FALSE), 
                    dze3 = Deriv(Z4[[1]],"z",cache.exp=FALSE), 
                    g1   = diffusion[1][[1]],
                    g3   = diffusion[3][[1]],  
                    C12=quote(C12),C13=quote(C13),C23=quote(C23),S1=quote(S1),S3=quote(S3)))))
    EqC1302 <- eval(Simplify(substitute(expression( 0.5 * ( g3 * dxxg1 + 2 * dxg1 * dxg3 + g1 * dxxg3)*S1 + 0.5 * ( g3 * dyyg1 + 2 * dyg1 * dyg3 + g1 * dyyg3)*S2 + 0.5 * ( g3 * dzzg1 + 2 * dzg1 * dzg3 + g1 * dzzg3)*S3 ), 
               list(g1   = diffusion[1][[1]],
                    g3   = diffusion[3][[1]], 
                    dxg1 = Deriv(diffusion[1],"x",cache.exp=FALSE)[[1]],
                    dxxg1= Deriv(diffusion[1],"x",nderiv=2,cache.exp=FALSE)[[1]],
                    dyg1 = Deriv(diffusion[1],"y",cache.exp=FALSE)[[1]],
                    dyyg1= Deriv(diffusion[1],"y",nderiv=2,cache.exp=FALSE)[[1]],
                    dzg1 = Deriv(diffusion[1],"z",cache.exp=FALSE)[[1]],
                    dzzg1= Deriv(diffusion[1],"z",nderiv=2,cache.exp=FALSE)[[1]],
                    dxg3 = Deriv(diffusion[3],"x",cache.exp=FALSE)[[1]],
                    dxxg3= Deriv(diffusion[3],"x",nderiv=2,cache.exp=FALSE)[[1]],
                    dyg3 = Deriv(diffusion[3],"y",cache.exp=FALSE)[[1]],
                    dyyg3= Deriv(diffusion[3],"y",nderiv=2,cache.exp=FALSE)[[1]],
                    dzg3 = Deriv(diffusion[3],"z",cache.exp=FALSE)[[1]],
                    dzzg3= Deriv(diffusion[3],"z",nderiv=2,cache.exp=FALSE)[[1]],
                    S1=quote(S1),S2=quote(S2),S3=quote(S3)))))
    EqC1303 <- eval(Simplify(substitute(expression( (g3*dxyg1+dyg1*dxg3+dxg1*dyg3+g1*dxyg3)*C12 + (g3*dxzg1+dzg1*dxg3+dxg1*dzg3+g1*dxzg3)*C13 + (g3*dyzg1+dzg1*dyg3+dyg1*dzg3+g1*dyzg3)*C23 ), 
               list(g1   = diffusion[1][[1]],
                    g3   = diffusion[3][[1]],
                    dxg1 = Deriv(diffusion[1],"x",cache.exp=FALSE)[[1]],
                    dyg1 = Deriv(diffusion[1],"y",cache.exp=FALSE)[[1]],
                    dzg1 = Deriv(diffusion[1],"z",cache.exp=FALSE)[[1]], 
                    dxg3 = Deriv(diffusion[3],"x",cache.exp=FALSE)[[1]],
                    dyg3 = Deriv(diffusion[3],"y",cache.exp=FALSE)[[1]],
                    dzg3 = Deriv(diffusion[3],"z",cache.exp=FALSE)[[1]],
                    dxyg1= Deriv(Deriv(diffusion[1],"x",cache.exp=FALSE),"y",cache.exp=FALSE)[[1]],
                    dxzg1= Deriv(Deriv(diffusion[1],"x",cache.exp=FALSE),"z",cache.exp=FALSE)[[1]],
                    dyzg1= Deriv(Deriv(diffusion[1],"y",cache.exp=FALSE),"z",cache.exp=FALSE)[[1]],
                    dxyg3= Deriv(Deriv(diffusion[3],"x",cache.exp=FALSE),"y",cache.exp=FALSE)[[1]],
                    dxzg3= Deriv(Deriv(diffusion[3],"x",cache.exp=FALSE),"z",cache.exp=FALSE)[[1]],
                    dyzg3= Deriv(Deriv(diffusion[3],"y",cache.exp=FALSE),"z",cache.exp=FALSE)[[1]],
                    C12=quote(C12),C13=quote(C13),C23=quote(C23)))))
    EqC1304 <- eval(Simplify(substitute(expression(e1 + e2 + e3 ), list(e1 = EqC1301[[1]], e2 = EqC1302[[1]],e3 = EqC1303[[1]]))))
    EqC130F <- eval(Simplify(substitute(substitute(e, list(x=quote(m1),y=quote(m2),z=quote(m3))), list(e = EqC1304[[1]]))))
    # C23(t)
    Z5 <- eval(Simplify(substitute(expression(e1 - 0.5 * (nu-1) * e2 * de2), list(e1 = drift[2][[1]], e2 = diffusion[2][[1]], de2 = Deriv(diffusion[2],"x",cache.exp=FALSE)[[1]],nu=nu))))
    Z6 <- eval(Simplify(substitute(expression(e1 - 0.5 * (nu-1) * e2 * de2), list(e1 = drift[3][[1]], e2 = diffusion[3][[1]], de2 = Deriv(diffusion[3],"z",cache.exp=FALSE)[[1]],nu=nu))))
    EqC2301 <- eval(Simplify(substitute(expression( dye3 * S2 + dze2 * S3 + ( dye2+dze3 ) *C23 + dxe3 * C12 + dxe2 * C13 + g2*g3), 
               list(dxe2 = Deriv(Z5[[1]],"x",cache.exp=FALSE),
                    dye2 = Deriv(Z5[[1]],"y",cache.exp=FALSE), 
                    dze2 = Deriv(Z5[[1]],"z",cache.exp=FALSE), 
                    dxe3 = Deriv(Z6[[1]],"x",cache.exp=FALSE), 
                    dye3 = Deriv(Z6[[1]],"y",cache.exp=FALSE), 
                    dze3 = Deriv(Z6[[1]],"z",cache.exp=FALSE), 
                    g2   = diffusion[2][[1]],
                    g3   = diffusion[3][[1]],  
                    C12=quote(C12),C13=quote(C13),C23=quote(C23),S2=quote(S2),S3=quote(S3)))))
    EqC2302 <- eval(Simplify(substitute(expression( 0.5 * ( g3 * dxxg2 + 2 * dxg2 * dxg3 + g2 * dxxg3)*S1 + 0.5 * ( g3 * dyyg2 + 2 * dyg2 * dyg3 + g2 * dyyg3)*S2 + 0.5 * ( g3 * dzzg2 + 2 * dzg2 * dzg3 + g2 * dzzg3)*S3 ), 
               list(g2   = diffusion[2][[1]],
                    g3   = diffusion[3][[1]], 
                    dxg2 = Deriv(diffusion[2],"x",cache.exp=FALSE)[[1]],
                    dxxg2= Deriv(diffusion[2],"x",nderiv=2,cache.exp=FALSE)[[1]],
                    dyg2 = Deriv(diffusion[2],"y",cache.exp=FALSE)[[1]],
                    dyyg2= Deriv(diffusion[2],"y",nderiv=2,cache.exp=FALSE)[[1]],
                    dzg2 = Deriv(diffusion[2],"z",cache.exp=FALSE)[[1]],
                    dzzg2= Deriv(diffusion[2],"z",nderiv=2,cache.exp=FALSE)[[1]],
                    dxg3 = Deriv(diffusion[3],"x",cache.exp=FALSE)[[1]],
                    dxxg3= Deriv(diffusion[3],"x",nderiv=2,cache.exp=FALSE)[[1]],
                    dyg3 = Deriv(diffusion[3],"y",cache.exp=FALSE)[[1]],
                    dyyg3= Deriv(diffusion[3],"y",nderiv=2,cache.exp=FALSE)[[1]],
                    dzg3 = Deriv(diffusion[3],"z",cache.exp=FALSE)[[1]],
                    dzzg3= Deriv(diffusion[3],"z",nderiv=2,cache.exp=FALSE)[[1]],
                    S1=quote(S1),S2=quote(S2),S3=quote(S3)))))
    EqC2303 <- eval(Simplify(substitute(expression( (g3*dxyg2+dyg2*dxg3+dxg2*dyg3+g2*dxyg3)*C12 + (g3*dxzg2+dzg2*dxg3+dxg2*dzg3+g2*dxzg3)*C13 + (g3*dyzg2+dzg2*dyg3+dyg2*dzg3+g2*dyzg3)*C23 ), 
               list(g2   = diffusion[2][[1]],
                    g3   = diffusion[3][[1]],
                    dxg2 = Deriv(diffusion[2],"x",cache.exp=FALSE)[[1]],
                    dyg2 = Deriv(diffusion[2],"y",cache.exp=FALSE)[[1]],
                    dzg2 = Deriv(diffusion[2],"z",cache.exp=FALSE)[[1]], 
                    dxg3 = Deriv(diffusion[3],"x",cache.exp=FALSE)[[1]],
                    dyg3 = Deriv(diffusion[3],"y",cache.exp=FALSE)[[1]],
                    dzg3 = Deriv(diffusion[3],"z",cache.exp=FALSE)[[1]],
                    dxyg2= Deriv(Deriv(diffusion[2],"x",cache.exp=FALSE),"y",cache.exp=FALSE)[[1]],
                    dxzg2= Deriv(Deriv(diffusion[2],"x",cache.exp=FALSE),"z",cache.exp=FALSE)[[1]],
                    dyzg2= Deriv(Deriv(diffusion[2],"y",cache.exp=FALSE),"z",cache.exp=FALSE)[[1]],
                    dxyg3= Deriv(Deriv(diffusion[3],"x",cache.exp=FALSE),"y",cache.exp=FALSE)[[1]],
                    dxzg3= Deriv(Deriv(diffusion[3],"x",cache.exp=FALSE),"z",cache.exp=FALSE)[[1]],
                    dyzg3= Deriv(Deriv(diffusion[3],"y",cache.exp=FALSE),"z",cache.exp=FALSE)[[1]],
                    C12=quote(C12),C13=quote(C13),C23=quote(C23)))))
    EqC2304 <- eval(Simplify(substitute(expression(e1 + e2 + e3 ), list(e1 = EqC2301[[1]], e2 = EqC2302[[1]],e3 = EqC2303[[1]]))))
    EqC230F <- eval(Simplify(substitute(substitute(e, list(x=quote(m1),y=quote(m2),z=quote(m3))), list(e = EqC2304[[1]]))))
    ## Solve ODE
    if (solve==TRUE){
    mod3d <- function(t,State,Pars){
       with(as.list(c(State, Pars)),{
          dm1  = eval(Eqm10F)
          dm2  = eval(Eqm20F)
          dm3  = eval(Eqm30F)
          dS1  = eval(EqS10F)
          dS2  = eval(EqS20F)
          dS3  = eval(EqS30F)
          dC12 = eval(EqC120F)
          dC13 = eval(EqC130F)
          dC23 = eval(EqC230F)
          return(list(c(dm1,dm2,dm3,dS1,dS2,dS3,dC12,dC13,dC23)))
          })
    }
    Sol <- deSolve::ode(func = mod3d, y = init,parms = parms, times = time,...)
	M = c(paste("| m1(",at,")  = ",sep=""),paste("| m2(",at,")  = ",sep=""),paste("| m3(",at,")  = ",sep=""))
    S = c(paste("| S1(",at,")  = ",sep=""),paste("| S2(",at,")  = ",sep=""),paste("| S3(",at,")  = ",sep=""))
    C = c(paste("| C12(",at,")  = ",sep=""),paste("| C13(",at,")  = ",sep=""),paste("| C23(",at,")  = ",sep=""))
	res1 = round(c(tail(Sol[,"m1"], n=1),tail(Sol[,"m2"], n=1),tail(Sol[,"m3"], n=1)),options()$digits)
    res2 = round(c(tail(Sol[,"S1"], n=1),tail(Sol[,"S2"], n=1),tail(Sol[,"S3"], n=1)),options()$digits)
    res3 = round(c(tail(Sol[,"C12"], n=1),tail(Sol[,"C13"], n=1),tail(Sol[,"C23"], n=1)),options()$digits)
	MEMF <- data.frame(M,res1,S,res2,C,res3,fix.empty.names = FALSE,row.names = c(" ", "","  "))}else{
    Sol <- NULL
	MEMF<- NULL}
    ## output
    structure(list(Means=list(m1=Eqm10F,m2=Eqm20F,m3=Eqm30F),Var=list(S1=EqS10F,S2=EqS20F,S3=EqS30F,C12=EqC120F,C13=EqC130F,C23=EqC230F),drift=drift,diffusion=diffusion,type=type,dim=length(drift),sol.ode=Sol,TAB=MEMF,call=match.call()),class="MEM.sde")
    }
}

print.MEM.sde <- function(x, digits=NULL, ...)
          {
    class(x) <- "MEM.sde"
	Ito = "It\xf4"
    Encoding(Ito) <- "latin1"
    if (!is.null(x$sol.ode)) {t0 = format(min(x$sol.ode[,"time"],na.rm = TRUE),digits=digits)
                              T  = format(max(x$sol.ode[,"time"],na.rm = TRUE),digits=digits)}else{
                              t0 = "t0";T = "T"
							  }
    if (x$dim==1){
      Dr <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)))), list(e = x$drift[[1]]))),width.cutoff=500L)   
      DD <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)))), list(e = x$diffusion[[1]]))),width.cutoff=500L)
      E  <- deparse(eval(substitute(substitute(e, list(m=quote(m(t)),S=quote(S(t)))), list(e = x$Means))),width.cutoff=500L)
      V  <- deparse(eval(substitute(substitute(e, list(m=quote(m(t)),S=quote(S(t)))), list(e = x$Var))),width.cutoff=500L)
      if(x$type=="ito"){
         cat(Ito," Sde 1D:","\n",        
             " | dX(t) = ", Dr," * dt + ", DD," * dW(t)","\n",
	         " | t in [",t0,",",T,"].","\n",      
             sep="")}else{
         cat("Stratonovich Sde 1D:","\n",
             " | dX(t) = ", Dr," * dt + ", DD," o dW(t)","\n",
	         " | t in [",t0,",",T,"].","\n", 
             sep="")}
    cat("\nMoment equations:","\n")
    cat(writeLines(strwrap(E, width = getOption("width"),indent = 0,exdent = 11, initial = " | dm(t) = ")))
    cat(writeLines(strwrap(V, width = getOption("width"),indent = 0,exdent = 11, initial = " | dS(t) = ")))
    if (!is.null(x$sol.ode)){
    cat("\nApproximation of moment at time ",tail(x$sol.ode[,"time"], n=1),"\n",
        " | m(",tail(x$sol.ode[,"time"], n=1),") = ", round(tail(x$sol.ode[,"m"], n=1),options()$digits) ,"\n", 
        " | S(",tail(x$sol.ode[,"time"], n=1),") = ", round(tail(x$sol.ode[,"S"], n=1),options()$digits) ,"\n",
    sep="")}
    }else if (x$dim==2){
	  Drx  <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)))), list(e = x$drift[1][[1]]))),width.cutoff=500L)   
      DDx  <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)))), list(e = x$diffusion[1][[1]]))),width.cutoff=500L)
	  Dry  <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)))), list(e = x$drift[2][[1]]))),width.cutoff=500L)   
      DDy  <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)))), list(e = x$diffusion[2][[1]]))),width.cutoff=500L)
      E1   <- deparse(eval(substitute(substitute(e, list(m1=quote(m1(t)),m2=quote(m2(t)),S1=quote(S1(t)),S2=quote(S2(t)),C12=quote(C12(t)))), list(e = x$Means$m1))),width.cutoff=500L)
      E2   <- deparse(eval(substitute(substitute(e, list(m1=quote(m1(t)),m2=quote(m2(t)),S1=quote(S1(t)),S2=quote(S2(t)),C12=quote(C12(t)))), list(e = x$Means$m2))),width.cutoff=500L)
      V1   <- deparse(eval(substitute(substitute(e, list(m1=quote(m1(t)),m2=quote(m2(t)),S1=quote(S1(t)),S2=quote(S2(t)),C12=quote(C12(t)))), list(e = x$Var$S1))),width.cutoff=500L)
      V2   <- deparse(eval(substitute(substitute(e, list(m1=quote(m1(t)),m2=quote(m2(t)),S1=quote(S1(t)),S2=quote(S2(t)),C12=quote(C12(t)))), list(e = x$Var$S2))),width.cutoff=500L)
      C12  <- deparse(eval(substitute(substitute(e, list(m1=quote(m1(t)),m2=quote(m2(t)),S1=quote(S1(t)),S2=quote(S2(t)),C12=quote(C12(t)))), list(e = x$Var$C12))),width.cutoff=500L)
      if(x$type=="ito"){
          cat(Ito," Sde 2D:","\n",
             " | dX(t) = ", Drx," * dt + ", DDx," * dW1(t)","\n", 
             " | dY(t) = ", Dry," * dt + ", DDy," * dW2(t)","\n",
	         " | t in [",t0,",",T,"].","\n",      
             sep="")}else{
          cat("Stratonovich Sde 2D:","\n",
             " | dX(t) = ", Drx," * dt + ", DDx," o dW1(t)","\n", 
             " | dY(t) = ", Dry," * dt + ", DDy," o dW2(t)","\n",
	         " | t in [",t0,",",T,"].","\n",      
             sep="")}
          cat("\nMoment equations:","\n")
          cat(writeLines(strwrap(E1, width = getOption("width"),indent = 0,exdent = 13, initial = " | dm1(t)  = ")), 
              writeLines(strwrap(E2, width = getOption("width"),indent = 0,exdent = 13, initial = " | dm2(t)  = ")), 
              writeLines(strwrap(V1, width = getOption("width"),indent = 0,exdent = 13, initial = " | dS1(t)  = ")),
              writeLines(strwrap(V2, width = getOption("width"),indent = 0,exdent = 13, initial = " | dS2(t)  = ")),
              writeLines(strwrap(C12,width = getOption("width"),indent = 0,exdent = 13, initial = " | dC12(t) = ")))
     if (!is.null(x$sol.ode)){
         cat("\nApproximation of moment at time ",tail(x$sol.ode[,"time"], n=1),sep="")
         print(x$TAB)           
                   }
    }else if (x$dim==3){
	  Drx  <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)),z=quote(Z(t)))), list(e = x$drift[1][[1]]))),width.cutoff=500L)   
      DDx  <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)),z=quote(Z(t)))), list(e = x$diffusion[1][[1]]))),width.cutoff=500L)
	  Dry  <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)),z=quote(Z(t)))), list(e = x$drift[2][[1]]))),width.cutoff=500L)   
      DDy  <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)),z=quote(Z(t)))), list(e = x$diffusion[2][[1]]))),width.cutoff=500L)
	  Drz  <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)),z=quote(Z(t)))), list(e = x$drift[3][[1]]))),width.cutoff=500L)   
      DDz  <- deparse(eval(substitute(substitute(e, list(x=quote(X(t)),y=quote(Y(t)),z=quote(Z(t)))), list(e = x$diffusion[3][[1]]))),width.cutoff=500L)
      E1   <- deparse(eval(substitute(substitute(e, list(m1=quote(m1(t)),m2=quote(m2(t)),m3=quote(m3(t)),S1=quote(S1(t)),S2=quote(S2(t)),S3=quote(S3(t)),C12=quote(C12(t)),C13=quote(C13(t)),C23=quote(C23(t)))), list(e = x$Means$m1))),width.cutoff=500L)
      E2   <- deparse(eval(substitute(substitute(e, list(m1=quote(m1(t)),m2=quote(m2(t)),m3=quote(m3(t)),S1=quote(S1(t)),S2=quote(S2(t)),S3=quote(S3(t)),C12=quote(C12(t)),C13=quote(C13(t)),C23=quote(C23(t)))), list(e = x$Means$m2))),width.cutoff=500L)
      E3   <- deparse(eval(substitute(substitute(e, list(m1=quote(m1(t)),m2=quote(m2(t)),m3=quote(m3(t)),S1=quote(S1(t)),S2=quote(S2(t)),S3=quote(S3(t)),C12=quote(C12(t)),C13=quote(C13(t)),C23=quote(C23(t)))), list(e = x$Means$m3))),width.cutoff=500L)	  
      V1   <- deparse(eval(substitute(substitute(e, list(m1=quote(m1(t)),m2=quote(m2(t)),m3=quote(m3(t)),S1=quote(S1(t)),S2=quote(S2(t)),S3=quote(S3(t)),C12=quote(C12(t)),C13=quote(C13(t)),C23=quote(C23(t)))), list(e = x$Var$S1))),width.cutoff=500L)
      V2   <- deparse(eval(substitute(substitute(e, list(m1=quote(m1(t)),m2=quote(m2(t)),m3=quote(m3(t)),S1=quote(S1(t)),S2=quote(S2(t)),S3=quote(S3(t)),C12=quote(C12(t)),C13=quote(C13(t)),C23=quote(C23(t)))), list(e = x$Var$S2))),width.cutoff=500L)
      V3   <- deparse(eval(substitute(substitute(e, list(m1=quote(m1(t)),m2=quote(m2(t)),m3=quote(m3(t)),S1=quote(S1(t)),S2=quote(S2(t)),S3=quote(S3(t)),C12=quote(C12(t)),C13=quote(C13(t)),C23=quote(C23(t)))), list(e = x$Var$S3))),width.cutoff=500L)
      C12  <- deparse(eval(substitute(substitute(e, list(m1=quote(m1(t)),m2=quote(m2(t)),m3=quote(m3(t)),S1=quote(S1(t)),S2=quote(S2(t)),S3=quote(S3(t)),C12=quote(C12(t)),C13=quote(C13(t)),C23=quote(C23(t)))), list(e = x$Var$C12))),width.cutoff=500L)
      C13  <- deparse(eval(substitute(substitute(e, list(m1=quote(m1(t)),m2=quote(m2(t)),m3=quote(m3(t)),S1=quote(S1(t)),S2=quote(S2(t)),S3=quote(S3(t)),C12=quote(C12(t)),C13=quote(C13(t)),C23=quote(C23(t)))), list(e = x$Var$C13))),width.cutoff=500L)
      C23  <- deparse(eval(substitute(substitute(e, list(m1=quote(m1(t)),m2=quote(m2(t)),m3=quote(m3(t)),S1=quote(S1(t)),S2=quote(S2(t)),S3=quote(S3(t)),C12=quote(C12(t)),C13=quote(C13(t)),C23=quote(C23(t)))), list(e = x$Var$C23))),width.cutoff=500L)
      if(x$type=="ito"){
         cat(Ito," Sde 3D:","\n",
             " | dX(t) = ", Drx," * dt + ", DDx," * dW1(t)","\n", 
             " | dY(t) = ", Dry," * dt + ", DDy," * dW2(t)","\n",
             " | dZ(t) = ", Drz," * dt + ", DDz," * dW3(t)","\n",
	         " | t in [",t0,",",T,"].","\n",     
             sep="")}else{
         cat("Stratonovich Sde 3D:","\n",
             " | dX(t) = ", Drx," * dt + ", DDx," o dW1(t)","\n", 
             " | dY(t) = ", Dry," * dt + ", DDy," o dW2(t)","\n",
             " | dZ(t) = ", Drz," * dt + ", DDz," o dW3(t)","\n",
	         " | t in [",t0,",",T,"].","\n",     
             sep="")} 
          cat("\nMoment equations:","\n")
          cat(writeLines(strwrap(E1, width = getOption("width"),indent = 0,exdent = 13, initial = " | dm1(t)  = ")), 
              writeLines(strwrap(E2, width = getOption("width"),indent = 0,exdent = 13, initial = " | dm2(t)  = ")),
              writeLines(strwrap(E3, width = getOption("width"),indent = 0,exdent = 13, initial = " | dm3(t)  = ")), 
              writeLines(strwrap(V1, width = getOption("width"),indent = 0,exdent = 13, initial = " | dS1(t)  = ")),
              writeLines(strwrap(V2, width = getOption("width"),indent = 0,exdent = 13, initial = " | dS2(t)  = ")),
              writeLines(strwrap(V3, width = getOption("width"),indent = 0,exdent = 13, initial = " | dS3(t)  = ")),
              writeLines(strwrap(C12, width = getOption("width"),indent = 0,exdent = 13, initial = " | dC12(t) = ")),
              writeLines(strwrap(C13, width = getOption("width"),indent = 0,exdent = 13, initial = " | dC13(t) = ")),
              writeLines(strwrap(C23, width = getOption("width"),indent = 0,exdent = 13, initial = " | dC23(t) = ")))
     if (!is.null(x$sol.ode)){
         cat("\nApproximation of moment at time ",tail(x$sol.ode[,"time"], n=1),sep="")
         print(x$TAB)
      }
    }
}

summary.MEM.sde <- function(object,at,...)
                {
    class(object) <- "MEM.sde"
    if (is.null(object$sol.ode)) {
	    stop("argument 'solve = FALSE'")
					   }else{
    if (missing(at)) {at = max(object$sol.ode[,"time"],na.rm = TRUE)}
    if (any(max(object$sol.ode[,"time"],na.rm = TRUE) < at || min(object$sol.ode[,"time"],na.rm = TRUE) > at) )  
	                    stop( " please use 'min(time) <= at <= max(time)'")
    #if (is.null(digits)){digits = base::options()$digits}
    if (object$dim==1){
       F   <- lapply(2:3,function(i) stats::approxfun(object$sol.ode[,"time"],object$sol.ode[,i]))
       Est <- sapply(1:length(F),function(i) F[[i]](at))
       cat("\nApproximation of moment at time ",at,"\n",
           " | m(",at,") = ", round(Est[1],options()$digits) ,"\n", 
           " | S(",at,") = ", round(Est[2],options()$digits) ,"\n",
           sep="")
    }else if (object$dim==2){
       F   <- lapply(2:6,function(i) stats::approxfun(object$sol.ode[,"time"],object$sol.ode[,i]))
       Est <- sapply(1:length(F),function(i) F[[i]](at))
	   M = c(paste("| m1(",at,")  = ",sep=""),paste("| m2(",at,")  = ",sep=""))
       S = c(paste("| S1(",at,")  = ",sep=""),paste("| S2(",at,")  = ",sep=""))
       C = c(paste("| C12(",at,")  = ",sep=""),"")
	   res1 = round(c(Est[1],Est[2]),options()$digits)
       res2 = round(c(Est[3],Est[4]),options()$digits)
       res3 = c(round(Est[5],options()$digits),"")
	   MEMF <- data.frame(M,res1,S,res2,C,res3,fix.empty.names = FALSE,row.names = c(" ", ""))
       cat("\nApproximation of moment at time ",at,sep="")
	   print(MEMF)
    } else if (object$dim==3){
       F   <- lapply(2:10,function(i) stats::approxfun(object$sol.ode[,"time"],object$sol.ode[,i]))
       Est <- sapply(1:length(F),function(i) F[[i]](at))
	   M = c(paste("| m1(",at,")  = ",sep=""),paste("| m2(",at,")  = ",sep=""),paste("| m3(",at,")  = ",sep=""))
       S = c(paste("| S1(",at,")  = ",sep=""),paste("| S2(",at,")  = ",sep=""),paste("| S3(",at,")  = ",sep=""))
       C = c(paste("| C12(",at,")  = ",sep=""),paste("| C13(",at,")  = ",sep=""),paste("| C23(",at,")  = ",sep=""))
	   res1 = round(c(Est[1],Est[2],Est[3]),options()$digits)
       res2 = round(c(Est[4],Est[5],Est[6]),options()$digits)
       res3 = round(c(Est[7],Est[8],Est[9]),options()$digits)
	   MEMF <- data.frame(M,res1,S,res2,C,res3,fix.empty.names = FALSE,row.names = c(" ", "","  "))
       cat("\nApproximation of moment at time ",at,sep="")
	   print(MEMF)
    } 
   }
}

