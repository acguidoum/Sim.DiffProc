## Fri Apr 03 08:54:53 2020
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



ABM <- function(N, ...)  UseMethod("ABM")

ABM.default <- function(N =1000,M=1,x0=0,t0=0,T=1,Dt=NULL,theta=1,sigma=1,...)
             {
			 
    if (!is.numeric(x0)) 
	   stop("'x0' must be numeric")
    if (any(!is.numeric(t0) || !is.numeric(T))) 
	   stop(" 't0' and 'T' must be numeric")
    if (any(!is.numeric(N)  || (N - floor(N) > 0) || N <= 1)) 
	   stop(" 'N' must be a positive integer ")
    if (any(!is.numeric(M)  || (M - floor(M) > 0) || M <= 0)) 
	   stop(" 'M' must be a positive integer ")
    if (any(t0 < 0 || T < 0 || T <= t0) ) 
	        stop(" please use positive times! (0 <= t0 < T) ")
    if (is.null(Dt)) {
        Dt <- (T - t0)/N
        t <- seq(t0, T, by=Dt)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
		T <- t[N + 1]
    }
    W <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
    X <- matrix(x0, N+1, M)
    for (i in 1L:N) {X[i + 1L,] <- X[i,] + theta * Dt + sigma * W[i,]  }
    name <- "X"
    name <- if(M > 1) paste("X",1:M,sep="")
    X <- ts(X, start = t0, deltat = Dt, names=name)
    return(X)
}

