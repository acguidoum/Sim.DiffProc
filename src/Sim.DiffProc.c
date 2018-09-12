/** 
 * Wed Aug 13 22:43:03 2014
 * The R Project for Statistical Computing
 * This file is part of the R package Sim.DiffProc
 * Original file Copyright © 2016 A.C. Guidoum <acguidoum@usthb.dz>; K. Boukhetala <kboukhetala@usthb.dz>
 * University of Science and Technology Houari Boumediene 
 * Department of Probabilities and Statistics
 * Faculty of Mathematics
 * BP 32 El-Alia, U.S.T.H.B, Algiers. Algeria.

 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * A copy of the GNU General Public License is available at <http://www.r-project.org/Licenses/>
 * Unlimited use and distribution (see LICENCE).

*/

/**
 * @file   Sim.DiffProc.c
 * @author A.C. Guidoum and K. Boukhetala
 * @date   2011-2014
 * @brief  C File for Numerical Methods for Multidimensional Stochastic Differential Equations. 
 */

#include "Sim.DiffProc.h"     			 
				   
/** Euler–Maruyama method (1-dim) 
 * @param  x0 initial value of the process at time t0
 * @param  t0 initial time
 * @param  delta time step of simulation (discretization)
 * @param  N number of simulation steps
 * @param  M number of trajectories
 * @funct  A drift coefficient
 * @funct  S diffusion coefficient
 * @param  rho the environtment on which to evaluate A and S
 * @return numerical solution (an time-series objects).
 * @author S.M. Iacus
 * @Ref S.M. Iacus: Simulation and Inference for Stochastic Differential Equations With R Examples, ISBN 978-0-387-75838-1, Springer, NY.
 */ 
 
SEXP Euler1d(SEXP x0, SEXP t0, SEXP delta, SEXP N, SEXP M,
             SEXP A, SEXP S, SEXP rho)
{
  SEXP X;
  double T, DELTA, sd, *rx0, *rX;
  double Z, itrx, a, s;
  int i, n, j, m;
  
  if(!isNumeric(x0))      error("`x0' must be numeric");
  if(!isNumeric(t0))      error("`t0' must be numeric");
  if(!isNumeric(delta))   error("`delta' must be numeric");
  if(!isInteger(N))       error("`N' must be integer");
  if(!isInteger(M))       error("`M' must be integer"); 
  if(!isFunction(A))      error("`A' must be a function");
  if(!isFunction(S))      error("`S' must be a function");
  if(!isEnvironment(rho)) error("`rho' must be an environment");
   
  PROTECT(x0 = AS_NUMERIC(x0));
  PROTECT(t0 = AS_NUMERIC(t0));
  PROTECT(delta = AS_NUMERIC(delta));
  PROTECT(N = AS_INTEGER(N));
  PROTECT(M = AS_INTEGER(M));
  n = *INTEGER(N);
  m = *INTEGER(M); 
  if(m>1)
   PROTECT(X = allocMatrix(REALSXP, n+1, m));
  else
   PROTECT(X = NEW_NUMERIC(n+1)); 
  rX = REAL(X);
  rx0 = REAL(x0);  
  for(j=0; j<m; j++)
   rX[j*(n+1)] = rx0[j]; 
   
  T = *REAL(t0);
  DELTA = *REAL(delta);
  sd = sqrt(DELTA);
  GetRNGstate();
    for(i=1; i< n+1; i++){
	 T = T + DELTA; 
     for(j=0; j< m; j++){
	  Z = rnorm(0,sd);
	  itrx = rX[i + j*(n+1) - 1]; 
	  a = fevalx(T,itrx,A,rho);
	  s = fevalx(T,itrx,S,rho);
      rX[i + (n+1)*j] = itrx + a*DELTA + s*Z;
	  }
	}
  PutRNGstate();
  UNPROTECT(6);
  return(X);
}

/** Euler–Maruyama method (2-dim)
 * @param  x0 initial value of the process X(t) at time t0
 * @param  y0 initial value of the process Y(t) at time t0
 * @param  t0 initial time
 * @param  delta time step of simulation (discretization)
 * @param  N number of simulation steps
 * @param  M number of trajectories
 * @funct  A1 drift coefficient of X(t)
 * @funct  S1 diffusion coefficient of X(t)
 * @funct  A2 drift coefficient of Y(t)
 * @funct  S2 diffusion coefficient of Y(t)
 * @param  rho the environtment on which to evaluate A1, S1, A2 and S2
 * @return numerical solution (an time-series objects).
 * @author A.C. Guidoum
 */ 
 
SEXP Euler2d(SEXP x0, SEXP y0, SEXP t0, SEXP delta, SEXP N, SEXP M,
             SEXP A1, SEXP S1, SEXP A2, SEXP S2, SEXP rho)
{
  SEXP X, Y, Val;
  double T, DELTA, sd, *rx0, *rX, *ry0, *rY, *rVal;
  double Z1, Z2, itrx, itry, a1, s1, a2, s2;
  int i, n, j, m;
  
  if(!isNumeric(x0))      error("`x0' must be numeric");
  if(!isNumeric(y0))      error("`y0' must be numeric");
  if(!isNumeric(t0))      error("`t0' must be numeric");
  if(!isNumeric(delta))   error("`delta' must be numeric");
  if(!isInteger(N))       error("`N' must be integer");
  if(!isInteger(M))       error("`M' must be integer"); 
  if(!isFunction(A1))     error("`A1' must be a function");
  if(!isFunction(S1))     error("`S1' must be a function");
  if(!isFunction(A2))     error("`A2' must be a function");
  if(!isFunction(S2))     error("`S2' must be a function");
  if(!isEnvironment(rho)) error("`rho' must be an environment");
   
  PROTECT(x0 = AS_NUMERIC(x0));
  PROTECT(y0 = AS_NUMERIC(y0));
  PROTECT(t0 = AS_NUMERIC(t0));
  PROTECT(delta = AS_NUMERIC(delta));
  PROTECT(N = AS_INTEGER(N));
  PROTECT(M = AS_INTEGER(M));
  n = *INTEGER(N);
  m = *INTEGER(M); 
  if(m>1){
   PROTECT(X = allocMatrix(REALSXP, n+1, m));
   PROTECT(Y = allocMatrix(REALSXP, n+1, m));
   PROTECT(Val = allocMatrix(REALSXP, n+1, 2*m));
        }
  else{
   PROTECT(X = NEW_NUMERIC(n+1)); 
   PROTECT(Y = NEW_NUMERIC(n+1)); 
   PROTECT(Val = allocMatrix(REALSXP, n+1, 2));
   }
  rX = REAL(X);
  rY = REAL(Y);
  rVal = REAL(Val);
  rx0 = REAL(x0); 
  ry0 = REAL(y0); 
  for(j=0; j<m; j++){
   rX[j*(n+1)] = rx0[j]; 
   rY[j*(n+1)] = ry0[j];
   rVal[(n+1)*j] = rx0[j];
   rVal[(n+1)*(j+m)] = ry0[j];
   }   
  T = *REAL(t0);
  DELTA = *REAL(delta);
  sd = sqrt(DELTA);
  GetRNGstate();
    for(i=1; i< n+1; i++){
	 T = T + DELTA; 
     for(j=0; j< m; j++){
	  Z1 = rnorm(0,sd);
	  Z2 = rnorm(0,sd);
	  itrx = rX[i + j*(n+1) - 1];
      itry = rY[i + j*(n+1) - 1]; 	  
	  a1 = fevalxy(T,itrx,itry,A1,rho);
	  s1 = fevalxy(T,itrx,itry,S1,rho);
	  a2 = fevalxy(T,itrx,itry,A2,rho);
	  s2 = fevalxy(T,itrx,itry,S2,rho);
      rX[i + (n+1)*j] = itrx + a1*DELTA + s1*Z1;
      rY[i + (n+1)*j] = itry + a2*DELTA + s2*Z2;
	  rVal[i+(n+1)*j] = rX[i+(n+1)*j];
	  rVal[i+(n+1)*(j+m)] = rY[i+(n+1)*j];
	  }
	}
  PutRNGstate();
  UNPROTECT(9);
  return Val;
}

/** Euler–Maruyama method (3-dim)
 * @param  x0 initial value of the process X(t) at time t0
 * @param  y0 initial value of the process Y(t) at time t0
 * @param  z0 initial value of the process Z(t) at time t0
 * @param  t0 initial time
 * @param  delta time step of simulation (discretization)
 * @param  N number of simulation steps
 * @param  M number of trajectories
 * @funct  A1 drift coefficient of X(t)
 * @funct  S1 diffusion coefficient of X(t)
 * @funct  A2 drift coefficient of Y(t)
 * @funct  S2 diffusion coefficient of Y(t)
 * @funct  A3 drift coefficient of Z(t)
 * @funct  S3 diffusion coefficient of Z(t)
 * @param  rho the environtment on which to evaluate A1, S1, A2, S2, A3 and S3
 * @return numerical solution (an time-series objects).
 * @author A.C. Guidoum
 */ 
 
SEXP Euler3d(SEXP x0, SEXP y0, SEXP z0, SEXP t0, SEXP delta, SEXP N, 
             SEXP M, SEXP A1, SEXP S1, SEXP A2, SEXP S2, SEXP A3, 
			 SEXP S3,SEXP rho)
{
  SEXP X, Y, Z, Val;
  double T, DELTA, sd, *rx0, *rX, *ry0, *rY, *rz0, *rZ, *rVal;
  double Z1, Z2, Z3, itrx, itry, itrz, a1, s1, a2, s2, a3, s3;
  int i, n, j, m;
  
  if(!isNumeric(x0))      error("`x0' must be numeric");
  if(!isNumeric(y0))      error("`y0' must be numeric");
  if(!isNumeric(z0))      error("`z0' must be numeric");
  if(!isNumeric(t0))      error("`t0' must be numeric");
  if(!isNumeric(delta))   error("`delta' must be numeric");
  if(!isInteger(N))       error("`N' must be integer");
  if(!isInteger(M))       error("`M' must be integer"); 
  if(!isFunction(A1))     error("`A1' must be a function");
  if(!isFunction(S1))     error("`S1' must be a function");
  if(!isFunction(A2))     error("`A2' must be a function");
  if(!isFunction(S2))     error("`S2' must be a function");
  if(!isFunction(A3))     error("`A3' must be a function");
  if(!isFunction(S3))     error("`S3' must be a function");
  if(!isEnvironment(rho)) error("`rho' must be an environment");
   
  PROTECT(x0 = AS_NUMERIC(x0));
  PROTECT(y0 = AS_NUMERIC(y0));
  PROTECT(z0 = AS_NUMERIC(z0));
  PROTECT(t0 = AS_NUMERIC(t0));
  PROTECT(delta = AS_NUMERIC(delta));
  PROTECT(N = AS_INTEGER(N));
  PROTECT(M = AS_INTEGER(M));
  n = *INTEGER(N);
  m = *INTEGER(M); 
  if(m>1){
   PROTECT(X = allocMatrix(REALSXP, n+1, m));
   PROTECT(Y = allocMatrix(REALSXP, n+1, m));
   PROTECT(Z = allocMatrix(REALSXP, n+1, m));
   PROTECT(Val = allocMatrix(REALSXP, n+1, 3*m));
        }
  else{
   PROTECT(X = NEW_NUMERIC(n+1)); 
   PROTECT(Y = NEW_NUMERIC(n+1)); 
   PROTECT(Z = NEW_NUMERIC(n+1)); 
   PROTECT(Val = allocMatrix(REALSXP, n+1, 3));
   }
  rX = REAL(X);
  rY = REAL(Y);
  rZ = REAL(Z);
  rVal = REAL(Val);
  rx0 = REAL(x0); 
  ry0 = REAL(y0); 
  rz0 = REAL(z0); 
  for(j=0; j<m; j++){
   rX[j*(n+1)] = rx0[j]; 
   rY[j*(n+1)] = ry0[j];
   rZ[j*(n+1)] = rz0[j];
   rVal[(n+1)*j] = rx0[j];
   rVal[(n+1)*(j+m)] = ry0[j];
   rVal[(n+1)*(j+2*m)] = rz0[j];
   }   
  T = *REAL(t0);
  DELTA = *REAL(delta);
  sd = sqrt(DELTA);
  GetRNGstate();
    for(i=1; i< n+1; i++){
	 T = T + DELTA; 
     for(j=0; j< m; j++){
	  Z1 = rnorm(0,sd);
	  Z2 = rnorm(0,sd);
	  Z3 = rnorm(0,sd);
	  itrx = rX[i + j*(n+1) - 1];
      itry = rY[i + j*(n+1) - 1]; 	
      itrz = rZ[i + j*(n+1) - 1]; 	  
	  a1 = fevalxyz(T,itrx,itry,itrz,A1,rho);
	  s1 = fevalxyz(T,itrx,itry,itrz,S1,rho);
	  a2 = fevalxyz(T,itrx,itry,itrz,A2,rho);
	  s2 = fevalxyz(T,itrx,itry,itrz,S2,rho);
	  a3 = fevalxyz(T,itrx,itry,itrz,A3,rho);
	  s3 = fevalxyz(T,itrx,itry,itrz,S3,rho);
      rX[i + (n+1)*j] = itrx + a1*DELTA + s1*Z1;
      rY[i + (n+1)*j] = itry + a2*DELTA + s2*Z2;
      rZ[i + (n+1)*j] = itrz + a3*DELTA + s3*Z3;
	  rVal[i+(n+1)*j] = rX[i+(n+1)*j];
	  rVal[i+(n+1)*(j+m)] = rY[i+(n+1)*j];
	  rVal[i+(n+1)*(j+2*m)] = rZ[i+(n+1)*j];
	  }
	}
  PutRNGstate();
  UNPROTECT(11);
  return Val;
}


/** Milstein method (1-dim)
 * @param  x0 initial value of the process at time t0
 * @param  t0 initial time
 * @param  delta time step of simulation (discretization)
 * @param  N number of simulation steps
 * @param  M number of trajectories
 * @funct  A drift coefficient
 * @funct  S diffusion coefficient
 * @funct  Sx partial derivative w.r.t. to x of S
 * @param  rho the environtment on which to evaluate A, S and Sx
 * @return numerical solution (an time-series objects).
 * @author S.M. Iacus
 * @Ref S.M. Iacus: Simulation and Inference for Stochastic Differential Equations With R Examples, ISBN 978-0-387-75838-1, Springer, NY. 
 */ 

SEXP Milstein1d(SEXP x0, SEXP t0, SEXP delta, SEXP N, SEXP M,
                SEXP A, SEXP S, SEXP Sx, SEXP rho)
{
  SEXP X;
  double T, DELTA, sd, *rx0, *rX;
  double Z, Z2, itrx, a, s, sx;
  int i, n, j, m;

  if(!isNumeric(x0))      error("`x0' must be numeric");
  if(!isNumeric(t0))      error("`t0' must be numeric");
  if(!isNumeric(delta))   error("`delta' must be numeric");
  if(!isInteger(N))       error("`N' must be integer");
  if(!isInteger(M))       error("`M' must be integer"); 
  if(!isFunction(A))      error("`A' must be a function");
  if(!isFunction(S))      error("`S' must be a function");
  if(!isFunction(Sx))      error("`Sx' must be a function");
  if(!isEnvironment(rho)) error("`rho' must be an environment");
   
  PROTECT(x0 = AS_NUMERIC(x0));
  PROTECT(t0 = AS_NUMERIC(t0));
  PROTECT(delta = AS_NUMERIC(delta));
  PROTECT(N = AS_INTEGER(N));
  PROTECT(M = AS_INTEGER(M));
  n = *INTEGER(N);
  m = *INTEGER(M); 
  if(m>1)
   PROTECT(X = allocMatrix(REALSXP, n+1, m));
  else
   PROTECT(X = NEW_NUMERIC(n+1)); 
  rX = REAL(X);
  rx0 = REAL(x0);  
  for(j=0; j<m; j++)
   rX[j*(n+1)] = rx0[j]; 
   
  T = *REAL(t0);
  DELTA = *REAL(delta);
  sd = sqrt(DELTA);
  GetRNGstate();
    for(i=1; i< n+1; i++){
	 T = T + DELTA; 
     for(j=0; j< m; j++){
	  Z = rnorm(0,sd);
	  Z2 = Z*Z;
	  itrx = rX[i + j*(n+1) - 1]; 
	  a = fevalx(T,itrx,A, rho);
	  s = fevalx(T,itrx,S, rho);
	  sx = fevalx(T,itrx,Sx, rho);
      rX[i + (n+1)*j] = itrx + a*DELTA + s*Z+ 0.5*s*sx*(Z2 - DELTA);						
	  }
	}
  PutRNGstate();
  UNPROTECT(6);
  return(X);
}

/** Milstein method (2-dim)
 * @param  x0 initial value of the process X(t) at time t0
 * @param  y0 initial value of the process Y(t) at time t0
 * @param  t0 initial time
 * @param  delta time step of simulation (discretization)
 * @param  N number of simulation steps
 * @param  M number of trajectories
 * @funct  A1, A2 drift coefficient
 * @funct  S1, S2 diffusion coefficient
 * @funct  S1x, S2x partial derivative w.r.t. to x (y) of S1 (S2)
 * @param  rho the environtment on which to evaluate A1, A2, S1, S2, S1x and S2x
 * @return numerical solution (an time-series objects).
 * @author A.C. Guidoum 
 */ 

SEXP Milstein2d(SEXP x0, SEXP y0, SEXP t0, SEXP delta, SEXP N, 
                SEXP M, SEXP A1, SEXP S1, SEXP S1x, SEXP A2, 
				SEXP S2, SEXP S2x, SEXP rho)
{
  SEXP X, Y, Val;
  double T, DELTA, sd, *rx0, *rX, *ry0, *rY, *rVal;
  double Z1, Z2, itrx, itry, a1, s1, a2, s2, s1x, s2x;
  int i, n, j, m;

  if(!isNumeric(x0))      error("`x0' must be numeric");
  if(!isNumeric(y0))      error("`y0' must be numeric");
  if(!isNumeric(t0))      error("`t0' must be numeric");
  if(!isNumeric(delta))   error("`delta' must be numeric");
  if(!isInteger(N))       error("`N' must be integer");
  if(!isInteger(M))       error("`M' must be integer"); 
  if(!isFunction(A1))     error("`A1' must be a function");
  if(!isFunction(S1))     error("`S1' must be a function");
  if(!isFunction(S1x))    error("`S1x' must be a function");
  if(!isFunction(A2))     error("`A2' must be a function");
  if(!isFunction(S2))     error("`S2' must be a function");
  if(!isFunction(S2x))    error("`S2x' must be a function");
  if(!isEnvironment(rho)) error("`rho' must be an environment");
   
  PROTECT(x0 = AS_NUMERIC(x0));
  PROTECT(y0 = AS_NUMERIC(y0));
  PROTECT(t0 = AS_NUMERIC(t0));
  PROTECT(delta = AS_NUMERIC(delta));
  PROTECT(N = AS_INTEGER(N));
  PROTECT(M = AS_INTEGER(M));
  n = *INTEGER(N);
  m = *INTEGER(M); 
    if(m>1){
   PROTECT(X = allocMatrix(REALSXP, n+1, m));
   PROTECT(Y = allocMatrix(REALSXP, n+1, m));
   PROTECT(Val = allocMatrix(REALSXP, n+1, 2*m));
        }
  else{
   PROTECT(X = NEW_NUMERIC(n+1)); 
   PROTECT(Y = NEW_NUMERIC(n+1)); 
   PROTECT(Val = allocMatrix(REALSXP, n+1, 2));
   }
  rX = REAL(X);
  rY = REAL(Y);
  rVal = REAL(Val);
  rx0 = REAL(x0); 
  ry0 = REAL(y0); 
  for(j=0; j<m; j++){
   rX[j*(n+1)] = rx0[j]; 
   rY[j*(n+1)] = ry0[j];
   rVal[(n+1)*j] = rx0[j];
   rVal[(n+1)*(j+m)] = ry0[j];
   }      
  T = *REAL(t0);
  DELTA = *REAL(delta);
  sd = sqrt(DELTA);
  GetRNGstate();
    for(i=1; i< n+1; i++){
	 T = T + DELTA; 
     for(j=0; j< m; j++){
	  Z1 = rnorm(0,sd);
	  Z2 = rnorm(0,sd);	  
	  itrx = rX[i + j*(n+1) - 1]; 
	  itry = rY[i + j*(n+1) - 1]; 
	  a1 = fevalxy(T,itrx,itry,A1, rho);
	  s1 = fevalxy(T,itrx,itry,S1, rho);
	  s1x = fevalxy(T,itrx,itry,S1x, rho);
	  a2 = fevalxy(T,itrx,itry,A2, rho);
	  s2 = fevalxy(T,itrx,itry,S2, rho);
	  s2x = fevalxy(T,itrx,itry,S2x, rho);
      rX[i + (n+1)*j] = itrx + a1*DELTA + s1*Z1+ 0.5*s1*s1x*(Z1*Z1 - DELTA);
      rY[i + (n+1)*j] = itry + a2*DELTA + s2*Z2+ 0.5*s2*s2x*(Z2*Z2 - DELTA);
	  rVal[i+(n+1)*j] = rX[i+(n+1)*j];
	  rVal[i+(n+1)*(j+m)] = rY[i+(n+1)*j];	  
	  }
	}
  PutRNGstate();
  UNPROTECT(9);
  return Val;
}

/** Milstein method (3-dim)
 * @param  x0 initial value of the process X(t) at time t0
 * @param  y0 initial value of the process Y(t) at time t0
 * @param  z0 initial value of the process Z(t) at time t0
 * @param  t0 initial time
 * @param  delta time step of simulation (discretization)
 * @param  N number of simulation steps
 * @param  M number of trajectories
 * @funct  A1, A2, A3 drift coefficient
 * @funct  S1, S2, S3 diffusion coefficient
 * @funct  S1x, S2x, S3x partial derivative w.r.t. to x (y)(z) of S1 (S2) (S3)
 * @param  rho the environtment on which to evaluate A1, A2, A3, S1, S2, S3, S1x, S2x and S3x
 * @return numerical solution (an time-series objects).
 * @author A.C. Guidoum 
 */ 

SEXP Milstein3d(SEXP x0, SEXP y0, SEXP z0, SEXP t0, SEXP delta, SEXP N, 
                SEXP M, SEXP A1, SEXP S1, SEXP S1x, SEXP A2, 
				SEXP S2, SEXP S2x, SEXP A3, SEXP S3, SEXP S3x, SEXP rho)
{
  SEXP X, Y, Z, Val;
  double T, DELTA, sd, *rx0, *rX, *ry0, *rY, *rz0, *rZ, *rVal;
  double Z1, Z2, Z3, itrx, itry, itrz, a1, s1, s1x, a2, s2, s2x, a3, s3, s3x;
  int i, n, j, m;

  if(!isNumeric(x0))      error("`x0' must be numeric");
  if(!isNumeric(y0))      error("`y0' must be numeric");
  if(!isNumeric(z0))      error("`z0' must be numeric");
  if(!isNumeric(t0))      error("`t0' must be numeric");
  if(!isNumeric(delta))   error("`delta' must be numeric");
  if(!isInteger(N))       error("`N' must be integer");
  if(!isInteger(M))       error("`M' must be integer"); 
  if(!isFunction(A1))     error("`A1' must be a function");
  if(!isFunction(S1))     error("`S1' must be a function");
  if(!isFunction(S1x))    error("`S1x' must be a function");
  if(!isFunction(A2))     error("`A2' must be a function");
  if(!isFunction(S2))     error("`S2' must be a function");
  if(!isFunction(S2x))    error("`S2x' must be a function");
  if(!isFunction(A3))     error("`A3' must be a function");
  if(!isFunction(S3))     error("`S3' must be a function");
  if(!isFunction(S3x))    error("`S3x' must be a function");
  if(!isEnvironment(rho)) error("`rho' must be an environment");
   
  PROTECT(x0 = AS_NUMERIC(x0));
  PROTECT(y0 = AS_NUMERIC(y0));
  PROTECT(z0 = AS_NUMERIC(z0));
  PROTECT(t0 = AS_NUMERIC(t0));
  PROTECT(delta = AS_NUMERIC(delta));
  PROTECT(N = AS_INTEGER(N));
  PROTECT(M = AS_INTEGER(M));
  n = *INTEGER(N);
  m = *INTEGER(M); 
    if(m>1){
   PROTECT(X = allocMatrix(REALSXP, n+1, m));
   PROTECT(Y = allocMatrix(REALSXP, n+1, m));
   PROTECT(Z = allocMatrix(REALSXP, n+1, m));
   PROTECT(Val = allocMatrix(REALSXP, n+1, 3*m));
        }
  else{
   PROTECT(X = NEW_NUMERIC(n+1)); 
   PROTECT(Y = NEW_NUMERIC(n+1)); 
   PROTECT(Z = NEW_NUMERIC(n+1));
   PROTECT(Val = allocMatrix(REALSXP, n+1, 3));
   }
  rX = REAL(X);
  rY = REAL(Y);
  rZ = REAL(Z);
  rVal = REAL(Val);
  rx0 = REAL(x0); 
  ry0 = REAL(y0); 
  rz0 = REAL(z0); 
  for(j=0; j<m; j++){
   rX[j*(n+1)] = rx0[j]; 
   rY[j*(n+1)] = ry0[j];
   rZ[j*(n+1)] = rz0[j];
   rVal[(n+1)*j] = rx0[j];
   rVal[(n+1)*(j+m)] = ry0[j];
   rVal[(n+1)*(j+2*m)] = rz0[j];
   }     
  T = *REAL(t0);
  DELTA = *REAL(delta);
  sd = sqrt(DELTA);
  GetRNGstate();
    for(i=1; i< n+1; i++){
	 T = T + DELTA; 
     for(j=0; j< m; j++){
	  Z1 = rnorm(0,sd);
	  Z2 = rnorm(0,sd);
	  Z3 = rnorm(0,sd);
	  itrx = rX[i + j*(n+1) - 1]; 
	  itry = rY[i + j*(n+1) - 1]; 
	  itrz = rZ[i + j*(n+1) - 1]; 
	  a1 = fevalxyz(T,itrx,itry,itrz,A1, rho);
	  s1 = fevalxyz(T,itrx,itry,itrz,S1, rho);
	  s1x = fevalxyz(T,itrx,itry,itrz,S1x, rho);
	  a2 = fevalxyz(T,itrx,itry,itrz,A2, rho);
	  s2 = fevalxyz(T,itrx,itry,itrz,S2, rho);
	  s2x = fevalxyz(T,itrx,itry,itrz,S2x, rho);
	  a3 = fevalxyz(T,itrx,itry,itrz,A3, rho);
	  s3 = fevalxyz(T,itrx,itry,itrz,S3, rho);
	  s3x = fevalxyz(T,itrx,itry,itrz,S3x, rho);
      rX[i + (n+1)*j] = itrx + a1*DELTA + s1*Z1+ 0.5*s1*s1x*(Z1*Z1 - DELTA);
      rY[i + (n+1)*j] = itry + a2*DELTA + s2*Z2+ 0.5*s2*s2x*(Z2*Z2 - DELTA);
	  rZ[i + (n+1)*j] = itrz + a3*DELTA + s3*Z3+ 0.5*s3*s3x*(Z3*Z3 - DELTA);
	  rVal[i+(n+1)*j] = rX[i+(n+1)*j];
	  rVal[i+(n+1)*(j+m)] = rY[i+(n+1)*j];
	  rVal[i+(n+1)*(j+2*m)] = rZ[i+(n+1)*j];	  
	  }
	}
  PutRNGstate();
  UNPROTECT(11);
  return Val;
}


/** Second Milstein method (1-dim)
 * @param  x0 initial value of the process at time t0
 * @param  t0 initial time
 * @param  delta time step of simulation (discretization)
 * @param  N number of simulation steps
 * @param  M number of trajectories
 * @funct  A drift coefficient
 * @funct  Ax partial derivative w.r.t. to x of A
 * @funct  Axx partial derivative w.r.t. to x of Ax
 * @funct  S diffusion coefficient
 * @funct  Sx partial derivative w.r.t. to x of S
 * @funct  Sxx partial derivative w.r.t. to x of Sx
 * @param  rho the environtment on which to evaluate A, Ax, Axx, S, Sx and Sxx
 * @return numerical solution (an time-series objects).
 * @author S.M. Iacus
 * @Ref S.M. Iacus: Simulation and Inference for Stochastic Differential Equations With R Examples, ISBN 978-0-387-75838-1, Springer, NY.
 */ 

SEXP SMilstein1d(SEXP x0, SEXP t0, SEXP delta, SEXP N, SEXP M,
                 SEXP A, SEXP Ax, SEXP Axx, SEXP S, SEXP Sx, 
				 SEXP Sxx, SEXP rho)
{
  SEXP X;
  double T, DELTA, sd, *rx0, *rX;
  double Z, Z2, itrx, a, ax, axx, s, sx , sxx;
  int i, n, j, m;

  if(!isNumeric(x0))      error("`x0' must be numeric");
  if(!isNumeric(t0))      error("`t0' must be numeric");
  if(!isNumeric(delta))   error("`delta' must be numeric");
  if(!isInteger(N))       error("`N' must be integer");
  if(!isInteger(M))       error("`M' must be integer"); 
  if(!isFunction(A))      error("`A' must be a function");
  if(!isFunction(Ax))     error("`Ax' must be a function");
  if(!isFunction(Axx))    error("`Axx' must be a function");
  if(!isFunction(S))      error("`S' must be a function");
  if(!isFunction(Sx))     error("`Sx' must be a function");
  if(!isFunction(Sxx))    error("`Sxx' must be a function");
  if(!isEnvironment(rho)) error("`rho' must be an environment");
   
  PROTECT(x0 = AS_NUMERIC(x0));
  PROTECT(t0 = AS_NUMERIC(t0));
  PROTECT(delta = AS_NUMERIC(delta));
  PROTECT(N = AS_INTEGER(N));
  PROTECT(M = AS_INTEGER(M));
  n = *INTEGER(N);
  m = *INTEGER(M); 
  if(m>1)
   PROTECT(X = allocMatrix(REALSXP, n+1, m));
  else
   PROTECT(X = NEW_NUMERIC(n+1)); 
  rX = REAL(X);
  rx0 = REAL(x0);  
  for(j=0; j<m; j++)
   rX[j*(n+1)] = rx0[j]; 
   
  T = *REAL(t0);
  DELTA = *REAL(delta);
  sd = sqrt(DELTA);
  GetRNGstate();
    for(i=1; i< n+1; i++){
	 T = T + DELTA; 
     for(j=0; j< m; j++){
	  Z = rnorm(0,sd);
	  Z2 = Z*Z;
	  itrx = rX[i + j*(n+1) - 1]; 
	  a = fevalx(T,itrx,A, rho);
	  ax = fevalx(T,itrx,Ax, rho);
	  axx = fevalx(T,itrx,Axx, rho);
	  s = fevalx(T,itrx,S, rho);
	  sx = fevalx(T,itrx,Sx, rho);
	  sxx = fevalx(T,itrx,Sxx, rho);
      rX[i + (n+1)*j] = itrx + a*DELTA + s*Z+ 0.5*s*sx*(Z2 - DELTA)+ pow(DELTA,1.5)*(0.5*a*sx+0.5*ax*s+0.25*s*s*sxx)*Z+
                        DELTA*DELTA* (0.5*a*ax+0.25*axx*s*s);	  
	  }
	}
  PutRNGstate();
  UNPROTECT(6);
  return(X);
}

/** Second Milstein method (2-dim)
 * @param  x0 initial value of the process X(t) at time t0
 * @param  y0 initial value of the process Y(t) at time t0
 * @param  t0 initial time
 * @param  delta time step of simulation (discretization)
 * @param  N number of simulation steps
 * @param  M number of trajectories
 * @funct  A1, A2 drift coefficient
 * @funct  A1x, A2x partial derivative w.r.t. to x (y) of A1 (A2)
 * @funct  A1xx, A2xx partial derivative w.r.t. to x (y) of A1x (A2x)
 * @funct  S1, S2 diffusion coefficient
 * @funct  S1x, S2x partial derivative w.r.t. to x (y) of S1 (S2)
 * @funct  S1xx, S2xx partial derivative w.r.t. to x (y) of S1x (S2x)
 * @param  rho the environtment on which to evaluate A1, A2, A1x, A2x, A1xx, A2xx, S1, S2, S1x, S2x, S1xx and S2xx
 * @return numerical solution (an time-series objects).
 * @author A.C. Guidoum 
 */ 

SEXP SMilstein2d(SEXP x0, SEXP y0, SEXP t0, SEXP delta, SEXP N, 
                 SEXP M, SEXP A1, SEXP A1x, SEXP A1xx, SEXP A2, 
				 SEXP A2x, SEXP A2xx, SEXP S1, SEXP S1x, SEXP S1xx, 
				 SEXP S2, SEXP S2x, SEXP S2xx, SEXP rho)
{
  SEXP X, Y, Val;
  double T, DELTA, sd, *rx0, *rX, *ry0, *rY, *rVal;
  double Z1, Z2, itrx, itry, a1, s1, a2, s2, a1x, a2x, a1xx, a2xx, s1x, s2x, s1xx, s2xx;
  int i, n, j, m;

  if(!isNumeric(x0))      error("`x0' must be numeric");
  if(!isNumeric(y0))      error("`y0' must be numeric");
  if(!isNumeric(t0))      error("`t0' must be numeric");
  if(!isNumeric(delta))   error("`delta' must be numeric");
  if(!isInteger(N))       error("`N' must be integer");
  if(!isInteger(M))       error("`M' must be integer"); 
  if(!isFunction(A1))     error("`A1' must be a function");
  if(!isFunction(A1x))    error("`A1x' must be a function");
  if(!isFunction(A1xx))   error("`A1xx' must be a function");
  if(!isFunction(S1))     error("`S1' must be a function");
  if(!isFunction(S1x))    error("`S1x' must be a function");
  if(!isFunction(S1xx))   error("`S1xx' must be a function");
  if(!isFunction(A2))     error("`A2' must be a function");
  if(!isFunction(A2x))    error("`A2x' must be a function");
  if(!isFunction(A2xx))   error("`A2xx' must be a function");
  if(!isFunction(S2))     error("`S2' must be a function");
  if(!isFunction(S2x))    error("`S2x' must be a function");
  if(!isFunction(S2xx))   error("`S2xx' must be a function");
  if(!isEnvironment(rho)) error("`rho' must be an environment");
   
  PROTECT(x0 = AS_NUMERIC(x0));
  PROTECT(y0 = AS_NUMERIC(y0));
  PROTECT(t0 = AS_NUMERIC(t0));
  PROTECT(delta = AS_NUMERIC(delta));
  PROTECT(N = AS_INTEGER(N));
  PROTECT(M = AS_INTEGER(M));
  n = *INTEGER(N);
  m = *INTEGER(M); 
    if(m>1){
   PROTECT(X = allocMatrix(REALSXP, n+1, m));
   PROTECT(Y = allocMatrix(REALSXP, n+1, m));
   PROTECT(Val = allocMatrix(REALSXP, n+1, 2*m));
        }
  else{
   PROTECT(X = NEW_NUMERIC(n+1)); 
   PROTECT(Y = NEW_NUMERIC(n+1)); 
   PROTECT(Val = allocMatrix(REALSXP, n+1, 2));
   }
  rX = REAL(X);
  rY = REAL(Y);
  rVal = REAL(Val);
  rx0 = REAL(x0); 
  ry0 = REAL(y0); 
  for(j=0; j<m; j++){
   rX[j*(n+1)] = rx0[j]; 
   rY[j*(n+1)] = ry0[j];
   rVal[(n+1)*j] = rx0[j];
   rVal[(n+1)*(j+m)] = ry0[j];
   }      
  T = *REAL(t0);
  DELTA = *REAL(delta);
  sd = sqrt(DELTA);
  GetRNGstate();
    for(i=1; i< n+1; i++){
	 T = T + DELTA; 
     for(j=0; j< m; j++){
	  Z1 = rnorm(0,sd);
	  Z2 = rnorm(0,sd);	  
	  itrx = rX[i + j*(n+1) - 1]; 
	  itry = rY[i + j*(n+1) - 1]; 
	  a1 = fevalxy(T,itrx,itry,A1, rho);
	  a1x = fevalxy(T,itrx,itry,A1x, rho);
	  a1xx = fevalxy(T,itrx,itry,A1xx, rho);
	  a2 = fevalxy(T,itrx,itry,A2, rho);
	  a2x = fevalxy(T,itrx,itry,A2x, rho);
	  a2xx = fevalxy(T,itrx,itry,A2xx, rho);
	  s1 = fevalxy(T,itrx,itry,S1, rho);
	  s1x = fevalxy(T,itrx,itry,S1x, rho);
	  s1xx = fevalxy(T,itrx,itry,S1xx, rho);
	  s2 = fevalxy(T,itrx,itry,S2, rho);
	  s2x = fevalxy(T,itrx,itry,S2x, rho);
	  s2xx = fevalxy(T,itrx,itry,S2xx, rho);
      rX[i + (n+1)*j] = itrx + a1*DELTA + s1*Z1+ 0.5*s1*s1x*(Z1*Z1 - DELTA)+ pow(DELTA,1.5)*(0.5*a1*s1x+0.5*a1x*s1+0.25*s1*s1*s1xx)*Z1+
                        DELTA*DELTA* (0.5*a1*a1x+0.25*a1xx*s1*s1);
      rY[i + (n+1)*j] = itry + a2*DELTA + s2*Z2+ 0.5*s2*s2x*(Z2*Z2 - DELTA)+ pow(DELTA,1.5)*(0.5*a2*s2x+0.5*a2x*s2+0.25*s2*s2*s2xx)*Z2+
                        DELTA*DELTA* (0.5*a2*a2x+0.25*a2xx*s2*s2);
	  rVal[i+(n+1)*j] = rX[i+(n+1)*j];
	  rVal[i+(n+1)*(j+m)] = rY[i+(n+1)*j];	  
	  }
	}
  PutRNGstate();
  UNPROTECT(9);
  return Val;
}

/** Second Milstein method (3-dim)
 * @param  x0 initial value of the process X(t) at time t0
 * @param  y0 initial value of the process Y(t) at time t0
 * @param  z0 initial value of the process Z(t) at time t0
 * @param  t0 initial time
 * @param  delta time step of simulation (discretization)
 * @param  N number of simulation steps
 * @param  M number of trajectories
 * @funct  A1, A2, A3 drift coefficient
 * @funct  A1x, A2x, A3x partial derivative w.r.t. to x (y)(z) of A1 (A2)(A3)
 * @funct  A1xx, A2xx, A3xx partial derivative w.r.t. to x (y)(z) of A1x (A2x)(A3x)
 * @funct  S1, S2, S3 diffusion coefficient
 * @funct  S1x, S2x, S3x partial derivative w.r.t. to x (y)(z) of S1 (S2)(S3)
 * @funct  S1xx, S2xx, S3xx partial derivative w.r.t. to x (y)(z) of S1x (S2x)(S3x)
 * @param  rho the environtment on which to evaluate A1, A2, A3, A1x, A2x, A3x, A1xx, A2xx, A3xx, S1, S2, S3, S1x, S2x, S3x, S1xx, S2xx and S3xx
 * @return numerical solution (an time-series objects).
 * @author A.C. Guidoum 
 */ 

SEXP SMilstein3d(SEXP x0, SEXP y0, SEXP z0, SEXP t0, SEXP delta, SEXP N, 
                 SEXP M, SEXP A1, SEXP A1x, SEXP A1xx, SEXP A2, 
				 SEXP A2x, SEXP A2xx, SEXP A3, SEXP A3x, SEXP A3xx, SEXP S1, 
				 SEXP S1x, SEXP S1xx, SEXP S2, SEXP S2x, SEXP S2xx, 
				 SEXP S3, SEXP S3x, SEXP S3xx, SEXP rho)
{
  SEXP X, Y, Z, Val;
  double T, DELTA, sd, *rx0, *rX, *ry0, *rY, *rz0, *rZ, *rVal;
  double Z1, Z2, Z3, itrx, itry, itrz, a1, s1, a2, s2, a3, s3, a1x, a2x, a3x, a1xx, a2xx, a3xx, s1x, s2x, s3x, s1xx, s2xx, s3xx;
  int i, n, j, m;

  if(!isNumeric(x0))      error("`x0' must be numeric");
  if(!isNumeric(y0))      error("`y0' must be numeric");
  if(!isNumeric(z0))      error("`z0' must be numeric");
  if(!isNumeric(t0))      error("`t0' must be numeric");
  if(!isNumeric(delta))   error("`delta' must be numeric");
  if(!isInteger(N))       error("`N' must be integer");
  if(!isInteger(M))       error("`M' must be integer"); 
  if(!isFunction(A1))     error("`A1' must be a function");
  if(!isFunction(A1x))    error("`A1x' must be a function");
  if(!isFunction(A1xx))   error("`A1xx' must be a function");
  if(!isFunction(S1))     error("`S1' must be a function");
  if(!isFunction(S1x))    error("`S1x' must be a function");
  if(!isFunction(S1xx))   error("`S1xx' must be a function");
  if(!isFunction(A2))     error("`A2' must be a function");
  if(!isFunction(A2x))    error("`A2x' must be a function");
  if(!isFunction(A2xx))   error("`A2xx' must be a function");
  if(!isFunction(S2))     error("`S2' must be a function");
  if(!isFunction(S2x))    error("`S2x' must be a function");
  if(!isFunction(S2xx))   error("`S2xx' must be a function");
  if(!isFunction(A3))     error("`A3' must be a function");
  if(!isFunction(A3x))    error("`A3x' must be a function");
  if(!isFunction(A3xx))   error("`A3xx' must be a function");
  if(!isFunction(S3))     error("`S3' must be a function");
  if(!isFunction(S3x))    error("`S3x' must be a function");
  if(!isFunction(S3xx))   error("`S3xx' must be a function");
  if(!isEnvironment(rho)) error("`rho' must be an environment");
   
  PROTECT(x0 = AS_NUMERIC(x0));
  PROTECT(y0 = AS_NUMERIC(y0));
  PROTECT(z0 = AS_NUMERIC(z0));
  PROTECT(t0 = AS_NUMERIC(t0));
  PROTECT(delta = AS_NUMERIC(delta));
  PROTECT(N = AS_INTEGER(N));
  PROTECT(M = AS_INTEGER(M));
  n = *INTEGER(N);
  m = *INTEGER(M); 
    if(m>1){
   PROTECT(X = allocMatrix(REALSXP, n+1, m));
   PROTECT(Y = allocMatrix(REALSXP, n+1, m));
   PROTECT(Z = allocMatrix(REALSXP, n+1, m));
   PROTECT(Val = allocMatrix(REALSXP, n+1, 3*m));
        }
  else{
   PROTECT(X = NEW_NUMERIC(n+1)); 
   PROTECT(Y = NEW_NUMERIC(n+1)); 
   PROTECT(Z = NEW_NUMERIC(n+1));
   PROTECT(Val = allocMatrix(REALSXP, n+1, 3));
   }
  rX = REAL(X);
  rY = REAL(Y);
  rZ = REAL(Z);
  rVal = REAL(Val);
  rx0 = REAL(x0); 
  ry0 = REAL(y0); 
  rz0 = REAL(z0); 
  for(j=0; j<m; j++){
   rX[j*(n+1)] = rx0[j]; 
   rY[j*(n+1)] = ry0[j];
   rZ[j*(n+1)] = rz0[j];
   rVal[(n+1)*j] = rx0[j];
   rVal[(n+1)*(j+m)] = ry0[j];
   rVal[(n+1)*(j+2*m)] = rz0[j];
   }      
  T = *REAL(t0);
  DELTA = *REAL(delta);
  sd = sqrt(DELTA);
  GetRNGstate();
    for(i=1; i< n+1; i++){
	 T = T + DELTA; 
     for(j=0; j< m; j++){
	  Z1 = rnorm(0,sd);
	  Z2 = rnorm(0,sd);
	  Z3 = rnorm(0,sd);
	  itrx = rX[i + j*(n+1) - 1]; 
	  itry = rY[i + j*(n+1) - 1]; 
	  itrz = rZ[i + j*(n+1) - 1]; 
	  a1 = fevalxyz(T,itrx,itry,itrz,A1, rho);
	  a1x = fevalxyz(T,itrx,itry,itrz,A1x, rho);
	  a1xx = fevalxyz(T,itrx,itry,itrz,A1xx, rho);
	  a2 = fevalxyz(T,itrx,itry,itrz,A2, rho);
	  a2x = fevalxyz(T,itrx,itry,itrz,A2x, rho);
	  a2xx = fevalxyz(T,itrx,itry,itrz,A2xx, rho);
	  a3 = fevalxyz(T,itrx,itry,itrz,A3, rho);
	  a3x = fevalxyz(T,itrx,itry,itrz,A3x, rho);
	  a3xx = fevalxyz(T,itrx,itry,itrz,A3xx, rho);
	  s1 = fevalxyz(T,itrx,itry,itrz,S1, rho);
	  s1x = fevalxyz(T,itrx,itry,itrz,S1x, rho);
	  s1xx = fevalxyz(T,itrx,itry,itrz,S1xx, rho);
	  s2 = fevalxyz(T,itrx,itry,itrz,S2, rho);
	  s2x = fevalxyz(T,itrx,itry,itrz,S2x, rho);
	  s2xx = fevalxyz(T,itrx,itry,itrz,S2xx, rho);
	  s3 = fevalxyz(T,itrx,itry,itrz,S3, rho);
	  s3x = fevalxyz(T,itrx,itry,itrz,S3x, rho);
	  s3xx = fevalxyz(T,itrx,itry,itrz,S3xx, rho);
      rX[i + (n+1)*j] = itrx + a1*DELTA + s1*Z1+ 0.5*s1*s1x*(Z1*Z1 - DELTA)+ pow(DELTA,1.5)*(0.5*a1*s1x+0.5*a1x*s1+0.25*s1*s1*s1xx)*Z1+
                        DELTA*DELTA* (0.5*a1*a1x+0.25*a1xx*s1*s1);
      rY[i + (n+1)*j] = itry + a2*DELTA + s2*Z2+ 0.5*s2*s2x*(Z2*Z2 - DELTA)+ pow(DELTA,1.5)*(0.5*a2*s2x+0.5*a2x*s2+0.25*s2*s2*s2xx)*Z2+
                        DELTA*DELTA* (0.5*a2*a2x+0.25*a2xx*s2*s2);
      rZ[i + (n+1)*j] = itrz + a3*DELTA + s3*Z3+ 0.5*s3*s3x*(Z3*Z3 - DELTA)+ pow(DELTA,1.5)*(0.5*a3*s3x+0.5*a3x*s3+0.25*s3*s3*s3xx)*Z3+
                        DELTA*DELTA* (0.5*a3*a3x+0.25*a3xx*s3*s3);
	  rVal[i+(n+1)*j] = rX[i+(n+1)*j];
	  rVal[i+(n+1)*(j+m)] = rY[i+(n+1)*j];
	  rVal[i+(n+1)*(j+2*m)] = rZ[i+(n+1)*j]; 
	  }
	}
  PutRNGstate();
  UNPROTECT(11);
  return Val;
}


/** Ito-Taylor Order 1.5 method (1-dim)
 * @param  x0 initial value of the process at time t0
 * @param  t0 initial time
 * @param  delta time step of simulation (discretization)
 * @param  N number of simulation steps
 * @param  M number of trajectories
 * @funct  A drift coefficient
 * @funct  Ax partial derivative w.r.t. to x of A
 * @funct  Axx partial derivative w.r.t. to x of Ax
 * @funct  S diffusion coefficient
 * @funct  Sx partial derivative w.r.t. to x of S
 * @funct  Sxx partial derivative w.r.t. to x of Sx
 * @param  rho the environtment on which to evaluate A, Ax, Axx, S, Sx and Sxx
 * @return numerical solution (an time-series objects).
 * @author S.M. Iacus
 * @Ref S.M. Iacus: Simulation and Inference for Stochastic Differential Equations With R Examples, ISBN 978-0-387-75838-1, Springer, NY.
 */ 

SEXP Sts1d(SEXP x0, SEXP t0, SEXP delta, SEXP N, SEXP M,
           SEXP A, SEXP Ax, SEXP Axx, SEXP S, SEXP Sx, 
		   SEXP Sxx, SEXP Z, SEXP U, SEXP rho)
{
  SEXP X;
  double T, DELTA, *rZ, *rU, *rx0, *rX;
  double z, u, itrx, a, ax, axx, s, sx , sxx;
  int i, n, j, m;

  if(!isNumeric(x0))      error("`x0' must be numeric");
  if(!isNumeric(t0))      error("`t0' must be numeric");
  if(!isNumeric(delta))   error("`delta' must be numeric");
  if(!isNumeric(Z))       error("`Z' must be numeric");
  if(!isNumeric(U))       error("`U' must be numeric");
  if(!isInteger(N))       error("`N' must be integer");
  if(!isInteger(M))       error("`M' must be integer"); 
  if(!isFunction(A))      error("`A' must be a function");
  if(!isFunction(Ax))     error("`Ax' must be a function");
  if(!isFunction(Axx))    error("`Axx' must be a function");
  if(!isFunction(S))      error("`S' must be a function");
  if(!isFunction(Sx))     error("`Sx' must be a function");
  if(!isFunction(Sxx))    error("`Sxx' must be a function");
  if(!isEnvironment(rho)) error("`rho' must be an environment");
   
  PROTECT(x0 = AS_NUMERIC(x0));
  PROTECT(t0 = AS_NUMERIC(t0));
  PROTECT(delta = AS_NUMERIC(delta));
  PROTECT(N = AS_INTEGER(N));
  PROTECT(M = AS_INTEGER(M));
  PROTECT(Z = AS_NUMERIC(Z));
  PROTECT(U = AS_NUMERIC(U));
  n = *INTEGER(N);
  m = *INTEGER(M); 
  if(m>1)
   PROTECT(X = allocMatrix(REALSXP, n+1, m));
  else
   PROTECT(X = NEW_NUMERIC(n+1)); 
  rX = REAL(X);
  rx0 = REAL(x0); 
  rZ = REAL(Z);
  rU = REAL(U);  
  for(j=0; j<m; j++)
   rX[j*(n+1)] = rx0[j]; 
   
  T = *REAL(t0);
  DELTA = *REAL(delta);
  GetRNGstate();
    for(i=1; i< n+1; i++){
	 T = T + DELTA; 
     for(j=0; j< m; j++){
	  itrx = rX[i + j*(n+1) - 1]; 
	  a = fevalx(T,itrx,A, rho);
	  ax = fevalx(T,itrx,Ax, rho);
	  axx = fevalx(T,itrx,Axx, rho);
	  s = fevalx(T,itrx,S, rho);
	  sx = fevalx(T,itrx,Sx, rho);
	  sxx = fevalx(T,itrx,Sxx, rho);	  
      z = rZ[j*n + i-1];
	  u = rU[j*n + i-1];
      rX[i + (n+1)*j] = itrx + a*DELTA + s*z+ 0.5*s*sx*(z*z - DELTA)+ ax*s*u+0.5*(a*ax+0.5*s*s*axx)*DELTA*DELTA+
						(a*sx+0.5*s*s*sxx)*(z*DELTA-u)+0.5*s*(s*sxx+sx*sx)*((z*z/3.0) -DELTA)*z;	  
	  }				
	}
  PutRNGstate();
  UNPROTECT(8);
  return(X);
}

/** Ito-Taylor Order 1.5 method (2-dim)
 * @param  x0 initial value of the process X(t) at time t0
 * @param  y0 initial value of the process Y(t) at time t0
 * @param  t0 initial time
 * @param  delta time step of simulation (discretization)
 * @param  N number of simulation steps
 * @param  M number of trajectories
 * @funct  A1, A2 drift coefficient
 * @funct  A1x, A2x partial derivative w.r.t. to x (y) of A1 (A2)
 * @funct  A1xx, A2xx partial derivative w.r.t. to x (y) of A1x (A2x)
 * @funct  S1, S2 diffusion coefficient
 * @funct  S1x, S2x partial derivative w.r.t. to x (y) of S1 (S2)
 * @funct  S1xx, S2xx partial derivative w.r.t. to x (y) of S1x (S2x)
 * @param  rho the environtment on which to evaluate A1, A2, A1x, A2x, A1xx, A2xx, S1, S2, S1x, S2x, S1xx and S2xx
 * @return numerical solution (an time-series objects).
 * @author A.C. Guidoum 
 */ 

SEXP Sts2d(SEXP x0, SEXP y0, SEXP t0, SEXP delta, SEXP N, SEXP M, 
           SEXP A1, SEXP A1x, SEXP A1xx, SEXP A2, SEXP A2x, SEXP A2xx, 
		   SEXP S1, SEXP S1x, SEXP S1xx, SEXP S2, SEXP S2x, SEXP S2xx, 
		   SEXP Z1, SEXP U1, SEXP Z2, SEXP U2, SEXP rho)
{
  SEXP X, Y, Val;
  double T, DELTA, *rZ1, *rU1, *rZ2, *rU2, *rx0, *rX, *ry0, *rY, *rVal;
  double z1, u1, z2, u2, itrx, itry, a1, s1, a2, s2, a1x, a2x, a1xx, a2xx, s1x, s2x, s1xx, s2xx;
  int i, n, j, m;

  if(!isNumeric(x0))      error("`x0' must be numeric");
  if(!isNumeric(y0))      error("`y0' must be numeric");
  if(!isNumeric(t0))      error("`t0' must be numeric");
  if(!isNumeric(delta))   error("`delta' must be numeric");
  if(!isNumeric(Z1))      error("`Z1' must be numeric");
  if(!isNumeric(U1))      error("`U1' must be numeric");
  if(!isNumeric(Z2))      error("`Z2' must be numeric");
  if(!isNumeric(U2))      error("`U2' must be numeric");
  if(!isInteger(N))       error("`N' must be integer");
  if(!isInteger(M))       error("`M' must be integer"); 
  if(!isFunction(A1))     error("`A1' must be a function");
  if(!isFunction(A1x))    error("`A1x' must be a function");
  if(!isFunction(A1xx))   error("`A1xx' must be a function");
  if(!isFunction(S1))     error("`S1' must be a function");
  if(!isFunction(S1x))    error("`S1x' must be a function");
  if(!isFunction(S1xx))   error("`S1xx' must be a function");
  if(!isFunction(A2))     error("`A2' must be a function");
  if(!isFunction(A2x))    error("`A2x' must be a function");
  if(!isFunction(A2xx))   error("`A2xx' must be a function");
  if(!isFunction(S2))     error("`S2' must be a function");
  if(!isFunction(S2x))    error("`S2x' must be a function");
  if(!isFunction(S2xx))   error("`S2xx' must be a function");
  if(!isEnvironment(rho)) error("`rho' must be an environment");
   
  PROTECT(x0 = AS_NUMERIC(x0));
  PROTECT(y0 = AS_NUMERIC(y0));
  PROTECT(t0 = AS_NUMERIC(t0));
  PROTECT(delta = AS_NUMERIC(delta));
  PROTECT(N = AS_INTEGER(N));
  PROTECT(M = AS_INTEGER(M));
  PROTECT(Z1 = AS_NUMERIC(Z1));
  PROTECT(U1 = AS_NUMERIC(U1));
  PROTECT(Z2 = AS_NUMERIC(Z2));
  PROTECT(U2 = AS_NUMERIC(U2));
  n = *INTEGER(N);
  m = *INTEGER(M); 
  if(m>1){
   PROTECT(X = allocMatrix(REALSXP, n+1, m));
   PROTECT(Y = allocMatrix(REALSXP, n+1, m));
   PROTECT(Val = allocMatrix(REALSXP, n+1, 2*m));
        }
  else{
   PROTECT(X = NEW_NUMERIC(n+1)); 
   PROTECT(Y = NEW_NUMERIC(n+1)); 
   PROTECT(Val = allocMatrix(REALSXP, n+1, 2));
   }
  rX = REAL(X);
  rY = REAL(Y);
  rZ1 = REAL(Z1);
  rU1 = REAL(U1); 
  rZ2 = REAL(Z2);
  rU2 = REAL(U2); 
  rVal = REAL(Val);
  rx0 = REAL(x0); 
  ry0 = REAL(y0); 
  for(j=0; j<m; j++){
   rX[j*(n+1)] = rx0[j]; 
   rY[j*(n+1)] = ry0[j];
   rVal[(n+1)*j] = rx0[j];
   rVal[(n+1)*(j+m)] = ry0[j];
   }      
  T = *REAL(t0);
  DELTA = *REAL(delta);
  GetRNGstate();
    for(i=1; i< n+1; i++){
	 T = T + DELTA; 
     for(j=0; j< m; j++){  
	  itrx = rX[i + j*(n+1) - 1]; 
	  itry = rY[i + j*(n+1) - 1]; 
	  a1 = fevalxy(T,itrx,itry,A1, rho);
	  a1x = fevalxy(T,itrx,itry,A1x, rho);
	  a1xx = fevalxy(T,itrx,itry,A1xx, rho);
	  a2 = fevalxy(T,itrx,itry,A2, rho);
	  a2x = fevalxy(T,itrx,itry,A2x, rho);
	  a2xx = fevalxy(T,itrx,itry,A2xx, rho);
	  s1 = fevalxy(T,itrx,itry,S1, rho);
	  s1x = fevalxy(T,itrx,itry,S1x, rho);
	  s1xx = fevalxy(T,itrx,itry,S1xx, rho);
	  s2 = fevalxy(T,itrx,itry,S2, rho);
	  s2x = fevalxy(T,itrx,itry,S2x, rho);
	  s2xx = fevalxy(T,itrx,itry,S2xx, rho);
      z1 = rZ1[j*n + i-1];
      z2 = rZ2[j*n + i-1];
	  u1 = rU1[j*n + i-1];
	  u2 = rU2[j*n + i-1];
      rX[i + (n+1)*j] = itrx + a1*DELTA + s1*z1+ 0.5*s1*s1x*(z1*z1 - DELTA)+ a1x*s1*u1+0.5*(a1*a1x+0.5*s1*s1*a1xx)*DELTA*DELTA+
						(a1*s1x+0.5*s1*s1*s1xx)*(z1*DELTA-u1)+0.5*s1*(s1*s1xx+s1x*s1x)*((z1*z1/3.0) -DELTA)*z1;
      rY[i + (n+1)*j] = itry + a2*DELTA + s2*z2+ 0.5*s2*s2x*(z2*z2 - DELTA)+ a2x*s2*u2+0.5*(a2*a2x+0.5*s2*s2*a2xx)*DELTA*DELTA+
						(a2*s2x+0.5*s2*s2*s2xx)*(z2*DELTA-u2)+0.5*s2*(s2*s2xx+s2x*s2x)*((z2*z2/3.0) -DELTA)*z2;
	  rVal[i+(n+1)*j] = rX[i+(n+1)*j];
	  rVal[i+(n+1)*(j+m)] = rY[i+(n+1)*j];	  
	  }
	}
  PutRNGstate();
  UNPROTECT(13);
  return Val;
}

/** Ito-Taylor Order 1.5 method (3-dim)
 * @param  x0 initial value of the process X(t) at time t0
 * @param  y0 initial value of the process Y(t) at time t0
 * @param  z0 initial value of the process Z(t) at time t0
 * @param  t0 initial time
 * @param  delta time step of simulation (discretization)
 * @param  N number of simulation steps
 * @param  M number of trajectories
 * @funct  A1, A2, A3 drift coefficient
 * @funct  A1x, A2x, A3x partial derivative w.r.t. to x (y)(z) of A1 (A2)(A3)
 * @funct  A1xx, A2xx, A3xx partial derivative w.r.t. to x (y)(z) of A1x (A2x)(A3x)
 * @funct  S1, S2, S3 diffusion coefficient
 * @funct  S1x, S2x, S3x partial derivative w.r.t. to x (y)(z) of S1 (S2)(S3)
 * @funct  S1xx, S2xx, S3xx partial derivative w.r.t. to x (y)(z) of S1x (S2x)(S3x)
 * @param  rho the environtment on which to evaluate A1, A2, A3, A1x, A2x, A3x, A1xx, A2xx, A3xx, S1, S2, S3, S1x, S2x, S3x, S1xx, S2xx and S3xx
 * @return numerical solution (an time-series objects).
 * @author A.C. Guidoum 
 */ 

SEXP Sts3d(SEXP x0, SEXP y0, SEXP z0, SEXP t0, SEXP delta, SEXP N, 
           SEXP M, SEXP A1, SEXP A1x, SEXP A1xx, SEXP A2, 
		   SEXP A2x, SEXP A2xx, SEXP A3, SEXP A3x, SEXP A3xx, SEXP S1, 
		   SEXP S1x, SEXP S1xx, SEXP S2, SEXP S2x, SEXP S2xx, 
		   SEXP S3, SEXP S3x, SEXP S3xx, SEXP Z1, SEXP U1, SEXP Z2, 
		   SEXP U2, SEXP Z3, SEXP U3, SEXP rho)
{
  SEXP X, Y, Z, Val;
  double T, DELTA, *rZ1, *rU1, *rZ2, *rU2, *rZ3, *rU3, *rx0, *rX, *ry0, *rY, *rz0, *rZ, *rVal;
  double z1, u1, z2, u2, z3, u3, itrx, itry, itrz, a1, s1, a2, s2, a3, s3, a1x, a2x, a3x, a1xx, a2xx, a3xx, s1x, s2x, s3x, s1xx, s2xx, s3xx;
  int i, n, j, m;

  if(!isNumeric(x0))      error("`x0' must be numeric");
  if(!isNumeric(y0))      error("`y0' must be numeric");
  if(!isNumeric(z0))      error("`z0' must be numeric");
  if(!isNumeric(t0))      error("`t0' must be numeric");
  if(!isNumeric(Z1))      error("`Z1' must be numeric");
  if(!isNumeric(U1))      error("`U1' must be numeric");
  if(!isNumeric(Z2))      error("`Z2' must be numeric");
  if(!isNumeric(U2))      error("`U2' must be numeric");
  if(!isNumeric(Z3))      error("`Z3' must be numeric");
  if(!isNumeric(U3))      error("`U3' must be numeric");
  if(!isNumeric(delta))   error("`delta' must be numeric");
  if(!isInteger(N))       error("`N' must be integer");
  if(!isInteger(M))       error("`M' must be integer"); 
  if(!isFunction(A1))     error("`A1' must be a function");
  if(!isFunction(A1x))    error("`A1x' must be a function");
  if(!isFunction(A1xx))   error("`A1xx' must be a function");
  if(!isFunction(S1))     error("`S1' must be a function");
  if(!isFunction(S1x))    error("`S1x' must be a function");
  if(!isFunction(S1xx))   error("`S1xx' must be a function");
  if(!isFunction(A2))     error("`A2' must be a function");
  if(!isFunction(A2x))    error("`A2x' must be a function");
  if(!isFunction(A2xx))   error("`A2xx' must be a function");
  if(!isFunction(S2))     error("`S2' must be a function");
  if(!isFunction(S2x))    error("`S2x' must be a function");
  if(!isFunction(S2xx))   error("`S2xx' must be a function");
  if(!isFunction(A3))     error("`A3' must be a function");
  if(!isFunction(A3x))    error("`A3x' must be a function");
  if(!isFunction(A3xx))   error("`A3xx' must be a function");
  if(!isFunction(S3))     error("`S3' must be a function");
  if(!isFunction(S3x))    error("`S3x' must be a function");
  if(!isFunction(S3xx))   error("`S3xx' must be a function");
  if(!isEnvironment(rho)) error("`rho' must be an environment");
   
  PROTECT(x0 = AS_NUMERIC(x0));
  PROTECT(y0 = AS_NUMERIC(y0));
  PROTECT(z0 = AS_NUMERIC(z0));
  PROTECT(t0 = AS_NUMERIC(t0));
  PROTECT(Z1 = AS_NUMERIC(Z1));
  PROTECT(U1 = AS_NUMERIC(U1));
  PROTECT(Z2 = AS_NUMERIC(Z2));
  PROTECT(U2 = AS_NUMERIC(U2));
  PROTECT(Z3 = AS_NUMERIC(Z3));
  PROTECT(U3 = AS_NUMERIC(U3));
  PROTECT(delta = AS_NUMERIC(delta));
  PROTECT(N = AS_INTEGER(N));
  PROTECT(M = AS_INTEGER(M));
  n = *INTEGER(N);
  m = *INTEGER(M); 
    if(m>1){
   PROTECT(X = allocMatrix(REALSXP, n+1, m));
   PROTECT(Y = allocMatrix(REALSXP, n+1, m));
   PROTECT(Z = allocMatrix(REALSXP, n+1, m));
   PROTECT(Val = allocMatrix(REALSXP, n+1, 3*m));
        }
  else{
   PROTECT(X = NEW_NUMERIC(n+1)); 
   PROTECT(Y = NEW_NUMERIC(n+1)); 
   PROTECT(Z = NEW_NUMERIC(n+1));
   PROTECT(Val = allocMatrix(REALSXP, n+1, 3));
   }
  rX = REAL(X);
  rY = REAL(Y);
  rZ = REAL(Z);
  rZ1 = REAL(Z1);
  rU1 = REAL(U1); 
  rZ2 = REAL(Z2);
  rU2 = REAL(U2); 
  rZ3 = REAL(Z3);
  rU3 = REAL(U3); 
  rVal = REAL(Val);
  rx0 = REAL(x0); 
  ry0 = REAL(y0); 
  rz0 = REAL(z0); 
  for(j=0; j<m; j++){
   rX[j*(n+1)] = rx0[j]; 
   rY[j*(n+1)] = ry0[j];
   rZ[j*(n+1)] = rz0[j];
   rVal[(n+1)*j] = rx0[j];
   rVal[(n+1)*(j+m)] = ry0[j];
   rVal[(n+1)*(j+2*m)] = rz0[j];
   }      
  T = *REAL(t0);
  DELTA = *REAL(delta);
  GetRNGstate();
    for(i=1; i< n+1; i++){
	 T = T + DELTA; 
     for(j=0; j< m; j++){  
	  itrx = rX[i + j*(n+1) - 1]; 
	  itry = rY[i + j*(n+1) - 1]; 
	  itrz = rZ[i + j*(n+1) - 1]; 
	  a1 = fevalxyz(T,itrx,itry,itrz,A1, rho);
	  a1x = fevalxyz(T,itrx,itry,itrz,A1x, rho);
	  a1xx = fevalxyz(T,itrx,itry,itrz,A1xx, rho);
	  a2 = fevalxyz(T,itrx,itry,itrz,A2, rho);
	  a2x = fevalxyz(T,itrx,itry,itrz,A2x, rho);
	  a2xx = fevalxyz(T,itrx,itry,itrz,A2xx, rho);
	  a3 = fevalxyz(T,itrx,itry,itrz,A3, rho);
	  a3x = fevalxyz(T,itrx,itry,itrz,A3x, rho);
	  a3xx = fevalxyz(T,itrx,itry,itrz,A3xx, rho);
	  s1 = fevalxyz(T,itrx,itry,itrz,S1, rho);
	  s1x = fevalxyz(T,itrx,itry,itrz,S1x, rho);
	  s1xx = fevalxyz(T,itrx,itry,itrz,S1xx, rho);
	  s2 = fevalxyz(T,itrx,itry,itrz,S2, rho);
	  s2x = fevalxyz(T,itrx,itry,itrz,S2x, rho);
	  s2xx = fevalxyz(T,itrx,itry,itrz,S2xx, rho);
	  s3 = fevalxyz(T,itrx,itry,itrz,S3, rho);
	  s3x = fevalxyz(T,itrx,itry,itrz,S3x, rho);
	  s3xx = fevalxyz(T,itrx,itry,itrz,S3xx, rho);
	  z1 = rZ1[j*n + i-1];
      z2 = rZ2[j*n + i-1];
      z3 = rZ3[j*n + i-1];
	  u1 = rU1[j*n + i-1];
	  u2 = rU2[j*n + i-1];
	  u3 = rU3[j*n + i-1];
      rX[i + (n+1)*j] = itrx + a1*DELTA + s1*z1+ 0.5*s1*s1x*(z1*z1 - DELTA)+ a1x*s1*u1+0.5*(a1*a1x+0.5*s1*s1*a1xx)*DELTA*DELTA+
						(a1*s1x+0.5*s1*s1*s1xx)*(z1*DELTA-u1)+0.5*s1*(s1*s1xx+s1x*s1x)*((z1*z1/3.0) -DELTA)*z1;
      rY[i + (n+1)*j] = itry + a2*DELTA + s2*z2+ 0.5*s2*s2x*(z2*z2 - DELTA)+ a2x*s2*u2+0.5*(a2*a2x+0.5*s2*s2*a2xx)*DELTA*DELTA+
						(a2*s2x+0.5*s2*s2*s2xx)*(z2*DELTA-u2)+0.5*s2*(s2*s2xx+s2x*s2x)*((z2*z2/3.0) -DELTA)*z2;
      rZ[i + (n+1)*j] = itrz + a3*DELTA + s3*z3+ 0.5*s3*s3x*(z3*z3 - DELTA)+ a3x*s3*u3+0.5*(a3*a3x+0.5*s3*s3*a3xx)*DELTA*DELTA+
						(a3*s3x+0.5*s3*s3*s3xx)*(z3*DELTA-u3)+0.5*s3*(s3*s3xx+s3x*s3x)*((z3*z3/3.0) -DELTA)*z3;
	  rVal[i+(n+1)*j] = rX[i+(n+1)*j];
	  rVal[i+(n+1)*(j+m)] = rY[i+(n+1)*j];
	  rVal[i+(n+1)*(j+2*m)] = rZ[i+(n+1)*j]; 
	  }
	}
  PutRNGstate();
  UNPROTECT(17);
  return Val;
}


/** Heun method (Order 2) (1-dim)
 * @param  x0 initial value of the process at time t0
 * @param  t0 initial time
 * @param  delta time step of simulation (discretization)
 * @param  N number of simulation steps
 * @param  M number of trajectories
 * @funct  A drift coefficient
 * @funct  S diffusion coefficient
 * @param  rho the environtment on which to evaluate A and S
 * @return numerical solution (an time-series objects).
 * @author A.C. Guidoum
 */ 
 
SEXP Heun1d(SEXP x0, SEXP t0, SEXP delta, SEXP N, SEXP M,
            SEXP A, SEXP S, SEXP rho)
{
  SEXP X, Y;
  double T, DELTA, sd, *rx0, *rX, *rY;
  double Z, itrx, ax, sx, ay, sy;
  int i, n, j, m;
  
  if(!isNumeric(x0))      error("`x0' must be numeric");
  if(!isNumeric(t0))      error("`t0' must be numeric");
  if(!isNumeric(delta))   error("`delta' must be numeric");
  if(!isInteger(N))       error("`N' must be integer");
  if(!isInteger(M))       error("`M' must be integer"); 
  if(!isFunction(A))      error("`A' must be a function");
  if(!isFunction(S))      error("`S' must be a function");
  if(!isEnvironment(rho)) error("`rho' must be an environment");
   
  PROTECT(x0 = AS_NUMERIC(x0));
  PROTECT(t0 = AS_NUMERIC(t0));
  PROTECT(delta = AS_NUMERIC(delta));
  PROTECT(N = AS_INTEGER(N));
  PROTECT(M = AS_INTEGER(M));
  n = *INTEGER(N);
  m = *INTEGER(M); 
  PROTECT(Y = NEW_NUMERIC(m));
  if(m>1)
   PROTECT(X = allocMatrix(REALSXP, n+1, m));
  else
   PROTECT(X = NEW_NUMERIC(n+1)); 
  rX = REAL(X);
  rY = REAL(Y);
  rx0 = REAL(x0);  
  for(j=0; j<m; j++)
   rX[j*(n+1)] = rx0[j]; 
   
  for(j=0; j<m; j++)
   rY[j] = rX[j*(n+1)];
   
  T = *REAL(t0);
  DELTA = *REAL(delta);
  sd = sqrt(DELTA);
  GetRNGstate();
    for(i=1; i< n+1; i++){
	 T = T + DELTA; 
     for(j=0; j< m; j++){	 
	  Z = rnorm(0,sd);
	  itrx = rX[i + j*(n+1) - 1]; 
	  ax = fevalx(T,itrx,A,rho);
	  sx = fevalx(T,itrx,S,rho);
	  rY[j] = itrx + ax*DELTA + sx*Z;
	  ay = fevalx(T,rY[j],A,rho);
	  sy = fevalx(T,rY[j],S,rho);
      rX[i + (n+1)*j] = itrx + 0.5 *(ax+ay)*DELTA + 0.5*(sx+sy)*Z;
	  }
	}
  PutRNGstate();
  UNPROTECT(7);
  return(X);
}

/** Heun method (Order 2) (2-dim)
 * @param  x0 initial value of the process X(t) at time t0
 * @param  y0 initial value of the process Y(t) at time t0
 * @param  t0 initial time
 * @param  delta time step of simulation (discretization)
 * @param  N number of simulation steps
 * @param  M number of trajectories
 * @funct  A1 drift coefficient of X(t)
 * @funct  S1 diffusion coefficient of X(t)
 * @funct  A2 drift coefficient of Y(t)
 * @funct  S2 diffusion coefficient of Y(t)
 * @param  rho the environtment on which to evaluate A1, S1, A2 and S2
 * @return numerical solution (an time-series objects).
 * @author A.C. Guidoum
 */ 
 
SEXP Heun2d(SEXP x0, SEXP y0, SEXP t0, SEXP delta, SEXP N, SEXP M,
            SEXP A1, SEXP S1, SEXP A2, SEXP S2, SEXP rho)
{
  SEXP X, Y, V1, V2, Val;
  double T, DELTA, sd, *rx0, *rX, *ry0, *rY, *rV1, *rV2, *rVal;
  double Z1, Z2, itrx, itry, a1x, s1x, a2x, s2x, a1v, s1v, a2v, s2v;
  int i, n, j, m;
  
  if(!isNumeric(x0))      error("`x0' must be numeric");
  if(!isNumeric(y0))      error("`y0' must be numeric");
  if(!isNumeric(t0))      error("`t0' must be numeric");
  if(!isNumeric(delta))   error("`delta' must be numeric");
  if(!isInteger(N))       error("`N' must be integer");
  if(!isInteger(M))       error("`M' must be integer"); 
  if(!isFunction(A1))     error("`A1' must be a function");
  if(!isFunction(S1))     error("`S1' must be a function");
  if(!isFunction(A2))     error("`A2' must be a function");
  if(!isFunction(S2))     error("`S2' must be a function");
  if(!isEnvironment(rho)) error("`rho' must be an environment");
   
  PROTECT(x0 = AS_NUMERIC(x0));
  PROTECT(y0 = AS_NUMERIC(y0));
  PROTECT(t0 = AS_NUMERIC(t0));
  PROTECT(delta = AS_NUMERIC(delta));
  PROTECT(N = AS_INTEGER(N));
  PROTECT(M = AS_INTEGER(M)); 
  n = *INTEGER(N);
  m = *INTEGER(M);
  PROTECT(V1 = NEW_NUMERIC(m));
  PROTECT(V2 = NEW_NUMERIC(m));   
  if(m>1){
   PROTECT(X = allocMatrix(REALSXP, n+1, m));
   PROTECT(Y = allocMatrix(REALSXP, n+1, m));
   PROTECT(Val = allocMatrix(REALSXP, n+1, 2*m));
        }
  else{
   PROTECT(X = NEW_NUMERIC(n+1)); 
   PROTECT(Y = NEW_NUMERIC(n+1)); 
   PROTECT(Val = allocMatrix(REALSXP, n+1, 2));
   }
  rX = REAL(X);
  rY = REAL(Y);
  rV1 = REAL(V1);
  rV2 = REAL(V2);
  rVal = REAL(Val);
  rx0 = REAL(x0); 
  ry0 = REAL(y0); 
  for(j=0; j<m; j++){
   rX[j*(n+1)] = rx0[j]; 
   rY[j*(n+1)] = ry0[j];
   rVal[(n+1)*j] = rx0[j];
   rVal[(n+1)*(j+m)] = ry0[j];
   }   
  for(j=0; j<m; j++){
   rV1[j] = rX[j*(n+1)];
   rV2[j] = rY[j*(n+1)];
   }
  T = *REAL(t0);
  DELTA = *REAL(delta);
  sd = sqrt(DELTA);
  GetRNGstate();
    for(i=1; i< n+1; i++){
	 T = T + DELTA; 
     for(j=0; j< m; j++){
	  Z1 = rnorm(0,sd);
	  Z2 = rnorm(0,sd); 
	  itrx = rX[i + j*(n+1) - 1];
      itry = rY[i + j*(n+1) - 1]; 
	  a1x = fevalxy(T,itrx,itry,A1,rho);
	  s1x = fevalxy(T,itrx,itry,S1,rho);
	  a2x = fevalxy(T,itrx,itry,A2,rho);
	  s2x = fevalxy(T,itrx,itry,S2,rho);
	  rV1[j] = itrx + a1x*DELTA + s1x*Z1;
	  rV2[j] = itry + a2x*DELTA + s2x*Z2;
	  a1v = fevalxy(T,rV1[j],itry,A1,rho);
	  s1v = fevalxy(T,rV1[j],itry,S1,rho);
	  a2v = fevalxy(T,itrx,rV2[j],A2,rho);
	  s2v = fevalxy(T,itrx,rV2[j],S2,rho);  
      rX[i + (n+1)*j] = itrx + 0.5 *(a1x+a1v)*DELTA + 0.5*(s1x+s1v)*Z1;
      rY[i + (n+1)*j] = itry + 0.5 *(a2x+a2v)*DELTA + 0.5*(s2x+s2v)*Z2;
	  rVal[i+(n+1)*j] = rX[i+(n+1)*j];
	  rVal[i+(n+1)*(j+m)] = rY[i+(n+1)*j];
	  }
	}
  PutRNGstate();
  UNPROTECT(11);
  return Val;
}

/** Heun method (Order 2) (3-dim)
 * @param  x0 initial value of the process X(t) at time t0
 * @param  y0 initial value of the process Y(t) at time t0
 * @param  z0 initial value of the process Z(t) at time t0
 * @param  t0 initial time
 * @param  delta time step of simulation (discretization)
 * @param  N number of simulation steps
 * @param  M number of trajectories
 * @funct  A1 drift coefficient of X(t)
 * @funct  S1 diffusion coefficient of X(t)
 * @funct  A2 drift coefficient of Y(t)
 * @funct  S2 diffusion coefficient of Y(t)
 * @funct  A3 drift coefficient of Z(t)
 * @funct  S3 diffusion coefficient of Z(t)
 * @param  rho the environtment on which to evaluate A1, S1, A2, S2, A3 and S3
 * @return numerical solution (an time-series objects).
 * @author A.C. Guidoum
 */ 
 
SEXP Heun3d(SEXP x0, SEXP y0, SEXP z0, SEXP t0, SEXP delta, SEXP N, 
            SEXP M, SEXP A1, SEXP S1, SEXP A2, SEXP S2, SEXP A3, 
			SEXP S3, SEXP rho)
{
  SEXP X, Y, Z, V1, V2, V3, Val;
  double T, DELTA, sd, *rx0, *rX, *ry0, *rY, *rz0, *rZ, *rV1, *rV2, *rV3, *rVal;
  double Z1, Z2, Z3, itrx, itry, itrz, a1x, s1x, a2x, s2x, a3x, s3x, a1v, s1v, a2v, s2v, a3v, s3v;
  int i, n, j, m;
  
  if(!isNumeric(x0))      error("`x0' must be numeric");
  if(!isNumeric(y0))      error("`y0' must be numeric");
  if(!isNumeric(z0))      error("`z0' must be numeric");
  if(!isNumeric(t0))      error("`t0' must be numeric");
  if(!isNumeric(delta))   error("`delta' must be numeric");
  if(!isInteger(N))       error("`N' must be integer");
  if(!isInteger(M))       error("`M' must be integer"); 
  if(!isFunction(A1))     error("`A1' must be a function");
  if(!isFunction(S1))     error("`S1' must be a function");
  if(!isFunction(A2))     error("`A2' must be a function");
  if(!isFunction(S2))     error("`S2' must be a function");
  if(!isFunction(A3))     error("`A3' must be a function");
  if(!isFunction(S3))     error("`S3' must be a function");
  if(!isEnvironment(rho)) error("`rho' must be an environment");
   
  PROTECT(x0 = AS_NUMERIC(x0));
  PROTECT(y0 = AS_NUMERIC(y0));
  PROTECT(z0 = AS_NUMERIC(z0));
  PROTECT(t0 = AS_NUMERIC(t0));
  PROTECT(delta = AS_NUMERIC(delta));
  PROTECT(N = AS_INTEGER(N));
  PROTECT(M = AS_INTEGER(M));
  n = *INTEGER(N);
  m = *INTEGER(M); 
  PROTECT(V1 = NEW_NUMERIC(m));
  PROTECT(V2 = NEW_NUMERIC(m)); 
  PROTECT(V3 = NEW_NUMERIC(m));
  if(m>1){
   PROTECT(X = allocMatrix(REALSXP, n+1, m));
   PROTECT(Y = allocMatrix(REALSXP, n+1, m));
   PROTECT(Z = allocMatrix(REALSXP, n+1, m));
   PROTECT(Val = allocMatrix(REALSXP, n+1, 3*m));
        }
  else{
   PROTECT(X = NEW_NUMERIC(n+1)); 
   PROTECT(Y = NEW_NUMERIC(n+1)); 
   PROTECT(Z = NEW_NUMERIC(n+1)); 
   PROTECT(Val = allocMatrix(REALSXP, n+1, 3));
   }
  rX = REAL(X);
  rY = REAL(Y);
  rZ = REAL(Z);
  rV1 = REAL(V1);
  rV2 = REAL(V2);
  rV3 = REAL(V3);
  rVal = REAL(Val);
  rx0 = REAL(x0); 
  ry0 = REAL(y0); 
  rz0 = REAL(z0); 
  for(j=0; j<m; j++){
   rX[j*(n+1)] = rx0[j]; 
   rY[j*(n+1)] = ry0[j];
   rZ[j*(n+1)] = rz0[j];
   rVal[(n+1)*j] = rx0[j];
   rVal[(n+1)*(j+m)] = ry0[j];
   rVal[(n+1)*(j+2*m)] = rz0[j];
   }   
  for(j=0; j<m; j++){
   rV1[j] = rX[j*(n+1)];
   rV2[j] = rY[j*(n+1)];
   rV3[j] = rZ[j*(n+1)];
   }
  T = *REAL(t0);
  DELTA = *REAL(delta);
  sd = sqrt(DELTA);
  GetRNGstate();
    for(i=1; i< n+1; i++){
	 T = T + DELTA; 
     for(j=0; j< m; j++){
	  Z1 = rnorm(0,sd);
	  Z2 = rnorm(0,sd); 
	  Z3 = rnorm(0,sd);
	  itrx = rX[i + j*(n+1) - 1];
      itry = rY[i + j*(n+1) - 1]; 	
      itrz = rZ[i + j*(n+1) - 1]; 	  
	  a1x = fevalxyz(T,itrx,itry,itrz,A1,rho);
	  s1x = fevalxyz(T,itrx,itry,itrz,S1,rho);
	  a2x = fevalxyz(T,itrx,itry,itrz,A2,rho);
	  s2x = fevalxyz(T,itrx,itry,itrz,S2,rho);
	  a3x = fevalxyz(T,itrx,itry,itrz,A3,rho);
	  s3x = fevalxyz(T,itrx,itry,itrz,S3,rho);
	  rV1[j] = itrx + a1x*DELTA + s1x*Z1;
	  rV2[j] = itry + a2x*DELTA + s2x*Z2;
	  rV3[j] = itrz + a3x*DELTA + s3x*Z3;
	  a1v = fevalxyz(T,rV1[j],itry,itrz,A1,rho);
	  s1v = fevalxyz(T,rV1[j],itry,itrz,S1,rho);
	  a2v = fevalxyz(T,itrx,rV2[j],itrz,A2,rho);
	  s2v = fevalxyz(T,itrx,rV2[j],itrz,S2,rho); 
	  a3v = fevalxyz(T,itrx,itry,rV3[j],A3,rho);
	  s3v = fevalxyz(T,itrx,itry,rV3[j],S3,rho);  	  
      rX[i + (n+1)*j] = itrx + 0.5 *(a1x+a1v)*DELTA + 0.5*(s1x+s1v)*Z1;
      rY[i + (n+1)*j] = itry + 0.5 *(a2x+a2v)*DELTA + 0.5*(s2x+s2v)*Z2;
      rZ[i + (n+1)*j] = itrz + 0.5 *(a3x+a3v)*DELTA + 0.5*(s3x+s3v)*Z3;
	  rVal[i+(n+1)*j] = rX[i+(n+1)*j];
	  rVal[i+(n+1)*(j+m)] = rY[i+(n+1)*j];
	  rVal[i+(n+1)*(j+2*m)] = rZ[i+(n+1)*j];
	  }
	}
  PutRNGstate();
  UNPROTECT(14);
  return Val;
}


/** Runge-Kutta method (Order 1, 2 and 3) (1-dim)
 * @param  x0 initial value of the process at time t0
 * @param  t0 initial time
 * @param  delta time step of simulation (discretization)
 * @param  N number of simulation steps
 * @param  M number of trajectories
 * @funct  A drift coefficient
 * @funct  S diffusion coefficient
 * @param  Order of convergence
 * @param  rho the environtment on which to evaluate A and S
 * @return numerical solution (an time-series objects).
 * @author A.C. Guidoum
 */ 
 
SEXP Rk1d(SEXP x0, SEXP t0, SEXP delta, SEXP N, SEXP M,
          SEXP A, SEXP S, SEXP Order, SEXP rho)
{
  SEXP X;
  double T, DELTA, sd, *rx0, *rX;
  double Z, itrx, a, a1, a2, s, s1, s2;
  int i, n, j, m, order;
  
  if(!isNumeric(x0))      error("`x0' must be numeric");
  if(!isNumeric(t0))      error("`t0' must be numeric");
  if(!isNumeric(delta))   error("`delta' must be numeric");
  if(!isInteger(N))       error("`N' must be integer");
  if(!isInteger(M))       error("`M' must be integer"); 
  if(!isInteger(Order))   error("`Order' must be integer"); 
  if(!isFunction(A))      error("`A' must be a function");
  if(!isFunction(S))      error("`S' must be a function");
  if(!isEnvironment(rho)) error("`rho' must be an environment");
   
  PROTECT(x0 = AS_NUMERIC(x0));
  PROTECT(t0 = AS_NUMERIC(t0));
  PROTECT(delta = AS_NUMERIC(delta));
  PROTECT(N = AS_INTEGER(N));
  PROTECT(M = AS_INTEGER(M));
  PROTECT(Order = AS_INTEGER(Order));
  n = *INTEGER(N);
  m = *INTEGER(M); 
  order = *INTEGER(Order);
  
  if(m>1)
   PROTECT(X = allocMatrix(REALSXP, n+1, m));
  else
   PROTECT(X = NEW_NUMERIC(n+1)); 
  rX = REAL(X);
  rx0 = REAL(x0);  
  for(j=0; j<m; j++)
   rX[j*(n+1)] = rx0[j]; 
   
  T = *REAL(t0);
  DELTA = *REAL(delta);
  sd = sqrt(DELTA);
  GetRNGstate();
  if (order==1){
    for(i=1; i< n+1; i++){
	 T = T + DELTA; 
     for(j=0; j< m; j++){	 
	  Z = rnorm(0,sd);
	  itrx = rX[i + j*(n+1) - 1]; 
	  a = fevalx(T,itrx,A,rho);
	  s = fevalx(T,itrx,S,rho);
	  s1 = fevalx(T+DELTA,itrx+s*sd,S,rho);
      rX[i + (n+1)*j] = itrx + a*DELTA + s*Z + (0.5/sd)*(s1-s)*(Z*Z-DELTA);	  
	  } 
	}
	}
   else if (order==2){
    for(i=1; i< n+1; i++){
	 T = T + DELTA; 
     for(j=0; j< m; j++){	 
	  Z = rnorm(0,sd);
	  itrx = rX[i + j*(n+1) - 1]; 
	  a = fevalx(T,itrx,A,rho);
	  s = fevalx(T,itrx,S,rho);
	  a1 = fevalx(T+DELTA,itrx+a*DELTA+s*Z,A,rho);	  
	  s1 = fevalx(T+DELTA,itrx+a*DELTA+sd*s,S,rho);
	  s2 = fevalx(T+DELTA,itrx+a*DELTA-sd*s,S,rho);  
      rX[i + (n+1)*j] = itrx + 0.5*(a+a1)*DELTA+0.25*(2*s+s1+s2)*Z+0.25*(s2-s1)*(sd - ((Z*Z)/sd));					
	  } 
	}
	}
   else if (order==3){
    for(i=1; i< n+1; i++){
	 T = T + DELTA; 
     for(j=0; j< m; j++){	 
	  Z = rnorm(0,sd);
	  itrx = rX[i + j*(n+1) - 1]; 
	  a = fevalx(T,itrx,A,rho);
	  s = fevalx(T,itrx,S,rho);
	  a1 = fevalx(T+0.5*DELTA,itrx+0.5*DELTA*a+s*Z,A,rho);	
	  s1 = fevalx(T+0.5*DELTA,itrx+0.5*DELTA*a+s*Z,S,rho);
	  a2 = fevalx(T+DELTA,itrx-a*DELTA+2*DELTA*a1+(2* s1-s)*Z,A,rho);	  	  
	  s2 = fevalx(T+DELTA,itrx-a*DELTA+2*DELTA*a1+(2* s1-s)*Z,S,rho); 
      rX[i + (n+1)*j] = itrx + (DELTA/6)*(a+4*a1+a2)+(Z/6)*(s+4*s1+s2);					
	  } 
	}
	}	
  PutRNGstate();
  UNPROTECT(7);
  return(X);
}

/** Runge-Kutta method (Order 1, 2 and 3) (2-dim)
 * @param  x0 initial value of the process X(t) at time t0
 * @param  y0 initial value of the process Y(t) at time t0
 * @param  t0 initial time
 * @param  delta time step of simulation (discretization)
 * @param  N number of simulation steps
 * @param  M number of trajectories
 * @funct  A1 drift coefficient of X(t)
 * @funct  S1 diffusion coefficient of X(t)
 * @funct  A2 drift coefficient of Y(t)
 * @funct  S2 diffusion coefficient of Y(t)
 * @param  Order of convergence
 * @param  rho the environtment on which to evaluate A1, S1, A2 and S2
 * @return numerical solution (an time-series objects).
 * @author A.C. Guidoum
 */ 
 
SEXP Rk2d(SEXP x0, SEXP y0, SEXP t0, SEXP delta, SEXP N, SEXP M,
          SEXP A1, SEXP S1, SEXP A2, SEXP S2, SEXP Order, SEXP rho)
{
  SEXP X, Y, Val;
  double T, DELTA, sd, *rx0, *rX, *ry0, *rY, *rVal;
  double Z1, Z2, itrx, itry, a1, a11, a12, s1, s11, s12, a2, a21, a22, s2, s21, s22;
  int i, n, j, m, order;
  
  if(!isNumeric(x0))      error("`x0' must be numeric");
  if(!isNumeric(y0))      error("`y0' must be numeric");
  if(!isNumeric(t0))      error("`t0' must be numeric");
  if(!isNumeric(delta))   error("`delta' must be numeric");
  if(!isInteger(N))       error("`N' must be integer");
  if(!isInteger(M))       error("`M' must be integer");
  if(!isInteger(Order))   error("`Order' must be integer");   
  if(!isFunction(A1))     error("`A1' must be a function");
  if(!isFunction(S1))     error("`S1' must be a function");
  if(!isFunction(A2))     error("`A2' must be a function");
  if(!isFunction(S2))     error("`S2' must be a function");
  if(!isEnvironment(rho)) error("`rho' must be an environment");
   
  PROTECT(x0 = AS_NUMERIC(x0));
  PROTECT(y0 = AS_NUMERIC(y0));
  PROTECT(t0 = AS_NUMERIC(t0));
  PROTECT(delta = AS_NUMERIC(delta));
  PROTECT(N = AS_INTEGER(N));
  PROTECT(M = AS_INTEGER(M)); 
  PROTECT(Order = AS_INTEGER(Order));
  n = *INTEGER(N);
  m = *INTEGER(M); 
  order = *INTEGER(Order); 
  if(m>1){
   PROTECT(X = allocMatrix(REALSXP, n+1, m));
   PROTECT(Y = allocMatrix(REALSXP, n+1, m));
   PROTECT(Val = allocMatrix(REALSXP, n+1, 2*m));
        }
  else{
   PROTECT(X = NEW_NUMERIC(n+1)); 
   PROTECT(Y = NEW_NUMERIC(n+1)); 
   PROTECT(Val = allocMatrix(REALSXP, n+1, 2));
   }
  rX = REAL(X);
  rY = REAL(Y);
  rVal = REAL(Val);
  rx0 = REAL(x0); 
  ry0 = REAL(y0); 
  for(j=0; j<m; j++){
   rX[j*(n+1)] = rx0[j]; 
   rY[j*(n+1)] = ry0[j];
   rVal[(n+1)*j] = rx0[j];
   rVal[(n+1)*(j+m)] = ry0[j];
   }   
  T = *REAL(t0);
  DELTA = *REAL(delta);
  sd = sqrt(DELTA);
  GetRNGstate();
   if (order==1){
    for(i=1; i< n+1; i++){
	 T = T + DELTA; 
     for(j=0; j< m; j++){
	  Z1 = rnorm(0,sd);
	  Z2 = rnorm(0,sd); 
	  itrx = rX[i + j*(n+1) - 1];
      itry = rY[i + j*(n+1) - 1]; 
	  a1 = fevalxy(T,itrx,itry,A1,rho);
	  s1 = fevalxy(T,itrx,itry,S1,rho);
	  s11 = fevalxy(T+DELTA,itrx+s1*sd,itry,S1,rho);
	  a2 = fevalxy(T,itrx,itry,A2,rho);
	  s2 = fevalxy(T,itrx,itry,S2,rho);
	  s21 = fevalxy(T+DELTA,itrx,itry+s2*sd,S2,rho);
      rX[i + (n+1)*j] = itrx + a1*DELTA + s1*Z1 + (0.5/sd)*(s11-s1)*(Z1*Z1-DELTA);
      rY[i + (n+1)*j] = itry + a2*DELTA + s2*Z2 + (0.5/sd)*(s21-s2)*(Z2*Z2-DELTA);	  
	  rVal[i+(n+1)*j] = rX[i+(n+1)*j];
	  rVal[i+(n+1)*(j+m)] = rY[i+(n+1)*j];
	  }
	 }
	}
   else if (order==2){
    for(i=1; i< n+1; i++){
	 T = T + DELTA; 
     for(j=0; j< m; j++){
	  Z1 = rnorm(0,sd);
	  Z2 = rnorm(0,sd); 
	  itrx = rX[i + j*(n+1) - 1];
      itry = rY[i + j*(n+1) - 1];   
	  a1 = fevalxy(T,itrx,itry,A1,rho);
	  s1 = fevalxy(T,itrx,itry,S1,rho);
	  a11 = fevalxy(T+DELTA,itrx+a1*DELTA+s1*Z1,itry,A1,rho);	  
	  s11 = fevalxy(T+DELTA,itrx+a1*DELTA+sd*s1,itry,S1,rho);
	  s12 = fevalxy(T+DELTA,itrx+a1*DELTA-sd*s1,itry,S1,rho);   
	  a2 = fevalxy(T,itrx,itry,A2,rho);
	  s2 = fevalxy(T,itrx,itry,S2,rho);
	  a21 = fevalxy(T+DELTA,itrx,itry+a2*DELTA+s2*Z2,A2,rho);	  
	  s21 = fevalxy(T+DELTA,itrx,itry+a2*DELTA+sd*s2,S2,rho);
	  s22 = fevalxy(T+DELTA,itrx,itry+a2*DELTA-sd*s2,S2,rho);   
      rX[i + (n+1)*j] = itrx + 0.5*(a1+a11)*DELTA+0.25*(2*s1+s11+s12)*Z1+0.25*(s12-s11)*(sd - ((Z1*Z1)/sd));	
      rY[i + (n+1)*j] = itry + 0.5*(a2+a21)*DELTA+0.25*(2*s2+s21+s22)*Z2+0.25*(s22-s21)*(sd - ((Z2*Z2)/sd));  
	  rVal[i+(n+1)*j] = rX[i+(n+1)*j];
	  rVal[i+(n+1)*(j+m)] = rY[i+(n+1)*j];
	  }
	 }
	}
   else if (order==3){
    for(i=1; i< n+1; i++){
	 T = T + DELTA; 
     for(j=0; j< m; j++){
	  Z1 = rnorm(0,sd);
	  Z2 = rnorm(0,sd); 
	  itrx = rX[i + j*(n+1) - 1];
      itry = rY[i + j*(n+1) - 1];   
	  a1 = fevalxy(T,itrx,itry,A1,rho);
	  s1 = fevalxy(T,itrx,itry,S1,rho);
	  a11 = fevalxy(T+0.5*DELTA,itrx+0.5*DELTA*a1+s1*Z1,itry,A1,rho);	
	  s11 = fevalxy(T+0.5*DELTA,itrx+0.5*DELTA*a1+s1*Z1,itry,S1,rho);
	  a12 = fevalxy(T+DELTA,itrx-a1*DELTA+2*DELTA*a11+(2* s11-s1)*Z1,itry,A1,rho);	  	  
	  s12 = fevalxy(T+DELTA,itrx-a1*DELTA+2*DELTA*a11+(2* s11-s1)*Z1,itry,S1,rho); 	  
	  a2 = fevalxy(T,itrx,itry,A2,rho);
	  s2 = fevalxy(T,itrx,itry,S2,rho);
	  a21 = fevalxy(T+0.5*DELTA,itrx,itry+0.5*DELTA*a2+s2*Z2,A2,rho);	
	  s21 = fevalxy(T+0.5*DELTA,itrx,itry+0.5*DELTA*a2+s2*Z2,S2,rho);
	  a22 = fevalxy(T+DELTA,itrx,itry-a2*DELTA+2*DELTA*a21+(2* s21-s2)*Z2,A2,rho);	  	  
	  s22 = fevalxy(T+DELTA,itrx,itry-a2*DELTA+2*DELTA*a21+(2* s21-s2)*Z2,S2,rho);   
      rX[i + (n+1)*j] = itrx + (DELTA/6)*(a1+4*a11+a12)+(Z1/6)*(s1+4*s11+s12);	
      rY[i + (n+1)*j] = itry + (DELTA/6)*(a2+4*a21+a22)+(Z2/6)*(s2+4*s21+s22);  
	  rVal[i+(n+1)*j] = rX[i+(n+1)*j];
	  rVal[i+(n+1)*(j+m)] = rY[i+(n+1)*j];
	  }
	 }
	}	
  PutRNGstate();
  UNPROTECT(10);
  return Val;
}

/** Runge-Kutta method (Order 1, 2 and 3) (3-dim)
 * @param  x0 initial value of the process X(t) at time t0
 * @param  y0 initial value of the process Y(t) at time t0
 * @param  z0 initial value of the process Z(t) at time t0
 * @param  t0 initial time
 * @param  delta time step of simulation (discretization)
 * @param  N number of simulation steps
 * @param  M number of trajectories
 * @funct  A1 drift coefficient of X(t)
 * @funct  S1 diffusion coefficient of X(t)
 * @funct  A2 drift coefficient of Y(t)
 * @funct  S2 diffusion coefficient of Y(t)
 * @funct  A3 drift coefficient of Z(t)
 * @funct  S3 diffusion coefficient of Z(t)
 * @param  rho the environtment on which to evaluate A1, S1, A2, S2, A3 and S3
 * @return numerical solution (an time-series objects).
 * @author A.C. Guidoum
 */ 
 
SEXP Rk3d(SEXP x0, SEXP y0, SEXP z0, SEXP t0, SEXP delta, SEXP N, 
            SEXP M, SEXP A1, SEXP S1, SEXP A2, SEXP S2, SEXP A3, 
			SEXP S3, SEXP Order, SEXP rho)
{
  SEXP X, Y, Z, Val;
  double T, DELTA, sd, *rx0, *rX, *ry0, *rY, *rz0, *rZ, *rVal;
  double Z1, Z2, Z3, itrx, itry, itrz, a1, a11, a12, s1, s11, s12, a2, a21, a22, s2, s21, s22, a3, a31, a32, s3, s31, s32;
  int i, n, j, m, order;
  
  if(!isNumeric(x0))      error("`x0' must be numeric");
  if(!isNumeric(y0))      error("`y0' must be numeric");
  if(!isNumeric(z0))      error("`z0' must be numeric");
  if(!isNumeric(t0))      error("`t0' must be numeric");
  if(!isNumeric(delta))   error("`delta' must be numeric");
  if(!isInteger(N))       error("`N' must be integer");
  if(!isInteger(M))       error("`M' must be integer"); 
  if(!isFunction(A1))     error("`A1' must be a function");
  if(!isFunction(S1))     error("`S1' must be a function");
  if(!isFunction(A2))     error("`A2' must be a function");
  if(!isFunction(S2))     error("`S2' must be a function");
  if(!isFunction(A3))     error("`A3' must be a function");
  if(!isFunction(S3))     error("`S3' must be a function");
  if(!isInteger(Order))   error("`Order' must be integer"); 
  if(!isEnvironment(rho)) error("`rho' must be an environment");
   
  PROTECT(x0 = AS_NUMERIC(x0));
  PROTECT(y0 = AS_NUMERIC(y0));
  PROTECT(z0 = AS_NUMERIC(z0));
  PROTECT(t0 = AS_NUMERIC(t0));
  PROTECT(delta = AS_NUMERIC(delta));
  PROTECT(N = AS_INTEGER(N));
  PROTECT(M = AS_INTEGER(M));
  PROTECT(Order = AS_INTEGER(Order));
  n = *INTEGER(N);
  m = *INTEGER(M); 
  order = *INTEGER(Order);
  if(m>1){
   PROTECT(X = allocMatrix(REALSXP, n+1, m));
   PROTECT(Y = allocMatrix(REALSXP, n+1, m));
   PROTECT(Z = allocMatrix(REALSXP, n+1, m));
   PROTECT(Val = allocMatrix(REALSXP, n+1, 3*m));
        }
  else{
   PROTECT(X = NEW_NUMERIC(n+1)); 
   PROTECT(Y = NEW_NUMERIC(n+1)); 
   PROTECT(Z = NEW_NUMERIC(n+1)); 
   PROTECT(Val = allocMatrix(REALSXP, n+1, 3));
   }
  rX = REAL(X);
  rY = REAL(Y);
  rZ = REAL(Z);
  rVal = REAL(Val);
  rx0 = REAL(x0); 
  ry0 = REAL(y0); 
  rz0 = REAL(z0); 
  for(j=0; j<m; j++){
   rX[j*(n+1)] = rx0[j]; 
   rY[j*(n+1)] = ry0[j];
   rZ[j*(n+1)] = rz0[j];
   rVal[(n+1)*j] = rx0[j];
   rVal[(n+1)*(j+m)] = ry0[j];
   rVal[(n+1)*(j+2*m)] = rz0[j];
   }   
  T = *REAL(t0);
  DELTA = *REAL(delta);
  sd = sqrt(DELTA);
  GetRNGstate();
   if (order==1){
    for(i=1; i< n+1; i++){
	 T = T + DELTA; 
     for(j=0; j< m; j++){
	  Z1 = rnorm(0,sd);
	  Z2 = rnorm(0,sd); 
	  Z3 = rnorm(0,sd);
	  itrx = rX[i + j*(n+1) - 1];
      itry = rY[i + j*(n+1) - 1]; 	
      itrz = rZ[i + j*(n+1) - 1]; 	
	  a1 = fevalxyz(T,itrx,itry,itrz,A1,rho);
	  s1 = fevalxyz(T,itrx,itry,itrz,S1,rho);
	  s11 = fevalxyz(T+DELTA,itrx+s1*sd,itry,itrz,S1,rho);
	  a2 = fevalxyz(T,itrx,itry,itrz,A2,rho);
	  s2 = fevalxyz(T,itrx,itry,itrz,S2,rho);
	  s21 = fevalxyz(T+DELTA,itrx,itry+s2*sd,itrz,S2,rho);  
	  a3 = fevalxyz(T,itrx,itry,itrz,A3,rho);
	  s3 = fevalxyz(T,itrx,itry,itrz,S3,rho);
	  s31 = fevalxyz(T+DELTA,itrx,itry,itrz+s3*sd,S3,rho);  
      rX[i + (n+1)*j] = itrx + a1*DELTA + s1*Z1 + (0.5/sd)*(s11-s1)*(Z1*Z1-DELTA);
      rY[i + (n+1)*j] = itry + a2*DELTA + s2*Z2 + (0.5/sd)*(s21-s2)*(Z2*Z2-DELTA);	 
      rZ[i + (n+1)*j] = itrz + a3*DELTA + s3*Z3 + (0.5/sd)*(s31-s3)*(Z3*Z3-DELTA);	  
	  rVal[i+(n+1)*j] = rX[i+(n+1)*j];
	  rVal[i+(n+1)*(j+m)] = rY[i+(n+1)*j];
	  rVal[i+(n+1)*(j+2*m)] = rZ[i+(n+1)*j];
	  }
	 }
	}
   else if (order==2){
    for(i=1; i< n+1; i++){
	 T = T + DELTA; 
     for(j=0; j< m; j++){
	  Z1 = rnorm(0,sd);
	  Z2 = rnorm(0,sd);
	  Z3 = rnorm(0,sd);
	  itrx = rX[i + j*(n+1) - 1];
      itry = rY[i + j*(n+1) - 1]; 	
      itrz = rZ[i + j*(n+1) - 1]; 	  
	  a1 = fevalxyz(T,itrx,itry,itrz,A1,rho);
	  s1 = fevalxyz(T,itrx,itry,itrz,S1,rho);
	  a11 = fevalxyz(T+DELTA,itrx+a1*DELTA+s1*Z1,itry,itrz,A1,rho);	  
	  s11 = fevalxyz(T+DELTA,itrx+a1*DELTA+sd*s1,itry,itrz,S1,rho);
	  s12 = fevalxyz(T+DELTA,itrx+a1*DELTA-sd*s1,itry,itrz,S1,rho);   
	  a2 = fevalxyz(T,itrx,itry,itrz,A2,rho);
	  s2 = fevalxyz(T,itrx,itry,itrz,S2,rho);
	  a21 = fevalxyz(T+DELTA,itrx,itry+a2*DELTA+s2*Z2,itrz,A2,rho);	  
	  s21 = fevalxyz(T+DELTA,itrx,itry+a2*DELTA+sd*s2,itrz,S2,rho);
	  s22 = fevalxyz(T+DELTA,itrx,itry+a2*DELTA-sd*s2,itrz,S2,rho); 
	  a3 = fevalxyz(T,itrx,itry,itrz,A3,rho);
	  s3 = fevalxyz(T,itrx,itry,itrz,S3,rho);
	  a31 = fevalxyz(T+DELTA,itrx,itry,itrz+a3*DELTA+s3*Z3,A3,rho);	  
	  s31 = fevalxyz(T+DELTA,itrx,itry,itrz+a3*DELTA+sd*s3,S3,rho);
	  s32 = fevalxyz(T+DELTA,itrx,itry,itrz+a3*DELTA-sd*s3,S3,rho);   	  
      rX[i + (n+1)*j] = itrx + 0.5*(a1+a11)*DELTA+0.25*(2*s1+s11+s12)*Z1+0.25*(s12-s11)*(sd - ((Z1*Z1)/sd));	
      rY[i + (n+1)*j] = itry + 0.5*(a2+a21)*DELTA+0.25*(2*s2+s21+s22)*Z2+0.25*(s22-s21)*(sd - ((Z2*Z2)/sd)); 
      rZ[i + (n+1)*j] = itrz + 0.5*(a3+a31)*DELTA+0.25*(2*s3+s31+s32)*Z3+0.25*(s32-s31)*(sd - ((Z3*Z3)/sd)); 	  
	  rVal[i+(n+1)*j] = rX[i+(n+1)*j];
	  rVal[i+(n+1)*(j+m)] = rY[i+(n+1)*j];
	  rVal[i+(n+1)*(j+2*m)] = rZ[i+(n+1)*j];
	  }
	 }
	}
   else if (order==3){
    for(i=1; i< n+1; i++){
	 T = T + DELTA; 
     for(j=0; j< m; j++){
	  Z1 = rnorm(0,sd);
	  Z2 = rnorm(0,sd);
	  Z3 = rnorm(0,sd);
	  itrx = rX[i + j*(n+1) - 1];
      itry = rY[i + j*(n+1) - 1]; 	
      itrz = rZ[i + j*(n+1) - 1]; 	   
	  a1 = fevalxyz(T,itrx,itry,itrz,A1,rho);
	  s1 = fevalxyz(T,itrx,itry,itrz,S1,rho);
	  a11 = fevalxyz(T+0.5*DELTA,itrx+0.5*DELTA*a1+s1*Z1,itry,itrz,A1,rho);	
	  s11 = fevalxyz(T+0.5*DELTA,itrx+0.5*DELTA*a1+s1*Z1,itry,itrz,S1,rho);
	  a12 = fevalxyz(T+DELTA,itrx-a1*DELTA+2*DELTA*a11+(2* s11-s1)*Z1,itry,itrz,A1,rho);	  	  
	  s12 = fevalxyz(T+DELTA,itrx-a1*DELTA+2*DELTA*a11+(2* s11-s1)*Z1,itry,itrz,S1,rho); 	  
	  a2 = fevalxyz(T,itrx,itry,itrz,A2,rho);
	  s2 = fevalxyz(T,itrx,itry,itrz,S2,rho);
	  a21 = fevalxyz(T+0.5*DELTA,itrx,itry+0.5*DELTA*a2+s2*Z2,itrz,A2,rho);	
	  s21 = fevalxyz(T+0.5*DELTA,itrx,itry+0.5*DELTA*a2+s2*Z2,itrz,S2,rho);
	  a22 = fevalxyz(T+DELTA,itrx,itry-a2*DELTA+2*DELTA*a21+(2* s21-s2)*Z2,itrz,A2,rho);	  	  
	  s22 = fevalxyz(T+DELTA,itrx,itry-a2*DELTA+2*DELTA*a21+(2* s21-s2)*Z2,itrz,S2,rho); 
	  a3 = fevalxyz(T,itrx,itry,itrz,A3,rho);
	  s3 = fevalxyz(T,itrx,itry,itrz,S3,rho);
	  a31 = fevalxyz(T+0.5*DELTA,itrx,itry,itrz+0.5*DELTA*a3+s3*Z3,A3,rho);	
	  s31 = fevalxyz(T+0.5*DELTA,itrx,itry,itrz+0.5*DELTA*a3+s3*Z3,S3,rho);
	  a32 = fevalxyz(T+DELTA,itrx,itry,itrz-a3*DELTA+2*DELTA*a31+(2* s31-s3)*Z3,A3,rho);	  	  
	  s32 = fevalxyz(T+DELTA,itrx,itry,itrz-a3*DELTA+2*DELTA*a31+(2* s31-s3)*Z3,S3,rho);	  
      rX[i + (n+1)*j] = itrx + (DELTA/6)*(a1+4*a11+a12)+(Z1/6)*(s1+4*s11+s12);	
      rY[i + (n+1)*j] = itry + (DELTA/6)*(a2+4*a21+a22)+(Z2/6)*(s2+4*s21+s22);  
      rZ[i + (n+1)*j] = itrz + (DELTA/6)*(a3+4*a31+a32)+(Z3/6)*(s3+4*s31+s32); 
	  rVal[i+(n+1)*j] = rX[i+(n+1)*j];
	  rVal[i+(n+1)*(j+m)] = rY[i+(n+1)*j];
	  rVal[i+(n+1)*(j+2*m)] = rZ[i+(n+1)*j];
	  }
	 }
	}
  PutRNGstate();
  UNPROTECT(12);
  return Val;
}


/** Predictor-Corrector method (1-dim)
 * @param  x0 initial value of the process at time t0
 * @param  t0 initial time
 * @param  delta time step of simulation (discretization)
 * @param  N number of simulation steps
 * @param  M number of trajectories
 * @param  alpha weight of the predictor-corrector scheme
 * @param  mu weight of the predictor-corrector scheme
 * @funct  A drift coefficient
 * @funct  S diffusion coefficient
 * @funct  Sx partial derivative w.r.t. to x of S
 * @param  rho the environtment on which to evaluate A, S and Sx
 * @return numerical solution (an time-series objects).
 * @author S.M. Iacus
 * @Ref S.M. Iacus: Simulation and Inference for Stochastic Differential Equations With R Examples, ISBN 978-0-387-75838-1, Springer, NY.
 */ 
 
SEXP Predcorr1d(SEXP x0, SEXP t0, SEXP delta, SEXP N, SEXP M,
                SEXP alpha, SEXP mu, SEXP A, SEXP S, SEXP Sx, SEXP rho)
{
  SEXP X;
  double T, DELTA, Alpha, Mu, sd, *rx0, *rX;
  double Z, itrx, a, s, s1, ss1, ss2;
  int i, n, j, m;
  
  if(!isNumeric(x0))      error("`x0' must be numeric");
  if(!isNumeric(t0))      error("`t0' must be numeric");
  if(!isNumeric(delta))   error("`delta' must be numeric");
  if(!isInteger(N))       error("`N' must be integer");
  if(!isInteger(M))       error("`M' must be integer"); 
  if(!isFunction(A))      error("`A' must be a function");
  if(!isFunction(S))      error("`S' must be a function");
  if(!isFunction(Sx))     error("`Sx' must be a function");
  if(!isNumeric(alpha))   error("`alpha' must be numeric");
  if(!isNumeric(mu))      error("`mu' must be numeric");
  if(!isEnvironment(rho)) error("`rho' must be an environment");
   
  PROTECT(x0 = AS_NUMERIC(x0));
  PROTECT(t0 = AS_NUMERIC(t0));
  PROTECT(delta = AS_NUMERIC(delta));
  PROTECT(N = AS_INTEGER(N));
  PROTECT(M = AS_INTEGER(M));
  PROTECT(alpha = AS_NUMERIC(alpha));
  PROTECT(mu = AS_NUMERIC(mu));

  n = *INTEGER(N);
  m = *INTEGER(M); 
  if(m>1)
   PROTECT(X = allocMatrix(REALSXP, n+1, m));
  else
   PROTECT(X = NEW_NUMERIC(n+1)); 
  rX = REAL(X);
  rx0 = REAL(x0);  
  Alpha = *REAL(alpha);
  Mu = *REAL(mu);
  for(j=0; j<m; j++)
   rX[j*(n+1)] = rx0[j]; 
   
  T = *REAL(t0);
  DELTA = *REAL(delta);
  sd = sqrt(DELTA);
  GetRNGstate();
    for(i=1; i< n+1; i++){
	 T = T + DELTA; 
     for(j=0; j< m; j++){
	  Z = rnorm(0,sd);
	  itrx = rX[i + j*(n+1) - 1]; 
	  a = fevalx(T,itrx,A,rho);
	  s = fevalx(T,itrx,S,rho);
	  s1 = fevalx(T+DELTA,itrx+a*DELTA+s*Z,S,rho);
	  ss1 = fevalx(T+DELTA,itrx+a*DELTA+s*Z,A,rho)- Mu *fevalx(T+DELTA,itrx+a*DELTA+s*Z,S,rho) * fevalx(T+DELTA,itrx+a*DELTA+s*Z,Sx,rho);
	  ss2 = fevalx(T,itrx,A,rho) - Mu *fevalx(T,itrx,S,rho) * fevalx(T,itrx,Sx,rho);
      rX[i + (n+1)*j] = itrx + (Alpha*ss1+(1-Alpha)*ss2)*DELTA + (Mu*s1+(1-Mu)*s)*Z;
	  }
	}
  PutRNGstate();
  UNPROTECT(8);
  return(X);
}
   
/** Predictor-Corrector method (2-dim)
 * @param  x0 initial value of the process X(t) at time t0
 * @param  y0 initial value of the process Y(t) at time t0
 * @param  t0 initial time
 * @param  delta time step of simulation (discretization)
 * @param  N number of simulation steps
 * @param  M number of trajectories
 * @param  alpha weight of the predictor-corrector scheme
 * @param  mu weight of the predictor-corrector scheme
 * @funct  A1, A2 drift coefficient
 * @funct  S1, S2 diffusion coefficient
 * @funct  S1x, S2x partial derivative w.r.t. to x (y) of S1 (S2)
 * @param  rho the environtment on which to evaluate A1, S1, S1x, A2, S2 and S2x
 * @return numerical solution (an time-series objects).
 * @author A.C. Guidoum
 */ 
 
SEXP Predcorr2d(SEXP x0, SEXP y0, SEXP t0, SEXP delta, SEXP N, SEXP M,
                SEXP alpha, SEXP mu, SEXP A1, SEXP S1, SEXP S1x,  
				SEXP A2, SEXP S2, SEXP S2x, SEXP rho)
{
  SEXP X, Y, Val;
  double T, DELTA, Alpha, Mu, sd, *rx0, *rX, *ry0, *rY, *rVal;
  double Z1, Z2, itrx, itry, a1, s1, s11, ss11, ss12, a2, s2, s21, ss21, ss22;
  int i, n, j, m;
  
  if(!isNumeric(x0))      error("`x0' must be numeric");
  if(!isNumeric(y0))      error("`y0' must be numeric");
  if(!isNumeric(t0))      error("`t0' must be numeric");
  if(!isNumeric(delta))   error("`delta' must be numeric");
  if(!isInteger(N))       error("`N' must be integer");
  if(!isInteger(M))       error("`M' must be integer"); 
  if(!isFunction(A1))     error("`A1' must be a function");
  if(!isFunction(S1))     error("`S1' must be a function");
  if(!isFunction(S1x))    error("`S1x' must be a function");
  if(!isFunction(A2))     error("`A2' must be a function");
  if(!isFunction(S2))     error("`S2' must be a function");
  if(!isFunction(S2x))    error("`S2x' must be a function");
  if(!isNumeric(alpha))   error("`alpha' must be numeric");
  if(!isNumeric(mu))      error("`mu' must be numeric");
  if(!isEnvironment(rho)) error("`rho' must be an environment");
   
  PROTECT(x0 = AS_NUMERIC(x0));
  PROTECT(y0 = AS_NUMERIC(y0));
  PROTECT(t0 = AS_NUMERIC(t0));
  PROTECT(delta = AS_NUMERIC(delta));
  PROTECT(N = AS_INTEGER(N));
  PROTECT(M = AS_INTEGER(M));
  PROTECT(alpha = AS_NUMERIC(alpha));
  PROTECT(mu = AS_NUMERIC(mu));

  n = *INTEGER(N);
  m = *INTEGER(M); 
  if(m>1){
   PROTECT(X = allocMatrix(REALSXP, n+1, m));
   PROTECT(Y = allocMatrix(REALSXP, n+1, m));
   PROTECT(Val = allocMatrix(REALSXP, n+1, 2*m));
        }
  else{
   PROTECT(X = NEW_NUMERIC(n+1)); 
   PROTECT(Y = NEW_NUMERIC(n+1)); 
   PROTECT(Val = allocMatrix(REALSXP, n+1, 2));
   }
  rX = REAL(X);
  rY = REAL(Y);
  rVal = REAL(Val);
  rx0 = REAL(x0); 
  ry0 = REAL(y0); 
  Alpha = *REAL(alpha);
  Mu = *REAL(mu);
  for(j=0; j<m; j++){
   rX[j*(n+1)] = rx0[j]; 
   rY[j*(n+1)] = ry0[j];
   rVal[(n+1)*j] = rx0[j];
   rVal[(n+1)*(j+m)] = ry0[j];
   }    
  T = *REAL(t0);
  DELTA = *REAL(delta);
  sd = sqrt(DELTA);
  GetRNGstate();
    for(i=1; i< n+1; i++){
	 T = T + DELTA; 
     for(j=0; j< m; j++){
	  Z1 = rnorm(0,sd);
	  Z2 = rnorm(0,sd); 
	  itrx = rX[i + j*(n+1) - 1];
      itry = rY[i + j*(n+1) - 1]; 	  
	  a1 = fevalxy(T,itrx,itry,A1,rho);
	  s1 = fevalxy(T,itrx,itry,S1,rho);
	  s11 = fevalxy(T+DELTA,itrx+a1*DELTA+s1*Z1,itry,S1,rho);
	  ss11 = fevalxy(T+DELTA,itrx+a1*DELTA+s1*Z1,itry,A1,rho)- Mu *fevalxy(T+DELTA,itrx+a1*DELTA+s1*Z1,itry,S1,rho) * fevalxy(T+DELTA,itrx+a1*DELTA+s1*Z1,itry,S1x,rho);
	  ss12 = fevalxy(T,itrx,itry,A1,rho) - Mu *fevalxy(T,itrx,itry,S1,rho) * fevalxy(T,itrx,itry,S1x,rho);
	  a2 = fevalxy(T,itrx,itry,A2,rho);
	  s2 = fevalxy(T,itrx,itry,S2,rho);
	  s21 = fevalxy(T+DELTA,itrx,itry+a2*DELTA+s2*Z2,S2,rho);
	  ss21 = fevalxy(T+DELTA,itrx,itry+a2*DELTA+s2*Z2,A2,rho)- Mu *fevalxy(T+DELTA,itrx,itry+a2*DELTA+s2*Z2,S2,rho) * fevalxy(T+DELTA,itrx,itry+a2*DELTA+s2*Z2,S2x,rho);
	  ss22 = fevalxy(T,itrx,itry,A2,rho) - Mu *fevalxy(T,itrx,itry,S2,rho) * fevalxy(T,itrx,itry,S2x,rho);
      rX[i + (n+1)*j] = itrx + (Alpha*ss11+(1-Alpha)*ss12)*DELTA + (Mu*s11+(1-Mu)*s1)*Z1;
      rY[i + (n+1)*j] = itry + (Alpha*ss21+(1-Alpha)*ss22)*DELTA + (Mu*s21+(1-Mu)*s2)*Z2;	  
	  rVal[i+(n+1)*j] = rX[i+(n+1)*j];
	  rVal[i+(n+1)*(j+m)] = rY[i+(n+1)*j];
	  }
	}
  PutRNGstate();
  UNPROTECT(11);
  return Val;
}
 
/** Predictor-Corrector method (3-dim)
 * @param  x0 initial value of the process X(t) at time t0
 * @param  y0 initial value of the process Y(t) at time t0
 * @param  z0 initial value of the process Z(t) at time t0
 * @param  t0 initial time
 * @param  delta time step of simulation (discretization)
 * @param  N number of simulation steps
 * @param  M number of trajectories
 * @funct  A1, A2, A3 drift coefficient
 * @funct  S1, S2, S3 diffusion coefficient
 * @funct  S1x, S2x, S3x partial derivative w.r.t. to x (y)(z) of S1 (S2) (S3)
 * @param  rho the environtment on which to evaluate A1, A2, A3, S1, S2, S3, S1x, S2x and S3x
 * @return numerical solution (an time-series objects).
 * @author A.C. Guidoum 
 */ 

SEXP Predcorr3d(SEXP x0, SEXP y0, SEXP z0, SEXP t0, SEXP delta, SEXP N, 
                SEXP M, SEXP alpha, SEXP mu, SEXP A1, SEXP S1, SEXP S1x, 
				SEXP A2, SEXP S2, SEXP S2x, SEXP A3, SEXP S3, SEXP S3x, 
				SEXP rho)
{
  SEXP X, Y, Z, Val;
  double T, DELTA, sd, Alpha, Mu, *rx0, *rX, *ry0, *rY, *rz0, *rZ, *rVal;
  double Z1, Z2, Z3, itrx, itry, itrz, a1, s1, s11, ss11, ss12, a2, s2, s21, ss21, ss22, a3, s3, s31, ss31, ss32;
  int i, n, j, m;

  if(!isNumeric(x0))      error("`x0' must be numeric");
  if(!isNumeric(y0))      error("`y0' must be numeric");
  if(!isNumeric(z0))      error("`z0' must be numeric");
  if(!isNumeric(t0))      error("`t0' must be numeric");
  if(!isNumeric(delta))   error("`delta' must be numeric");
  if(!isInteger(N))       error("`N' must be integer");
  if(!isInteger(M))       error("`M' must be integer"); 
  if(!isNumeric(alpha))   error("`alpha' must be numeric");
  if(!isNumeric(mu))      error("`mu' must be numeric");
  if(!isFunction(A1))     error("`A1' must be a function");
  if(!isFunction(S1))     error("`S1' must be a function");
  if(!isFunction(S1x))    error("`S1x' must be a function");
  if(!isFunction(A2))     error("`A2' must be a function");
  if(!isFunction(S2))     error("`S2' must be a function");
  if(!isFunction(S2x))    error("`S2x' must be a function");
  if(!isFunction(A3))     error("`A3' must be a function");
  if(!isFunction(S3))     error("`S3' must be a function");
  if(!isFunction(S3x))    error("`S3x' must be a function");
  if(!isEnvironment(rho)) error("`rho' must be an environment");
   
  PROTECT(x0 = AS_NUMERIC(x0));
  PROTECT(y0 = AS_NUMERIC(y0));
  PROTECT(z0 = AS_NUMERIC(z0));
  PROTECT(t0 = AS_NUMERIC(t0));
  PROTECT(delta = AS_NUMERIC(delta));
  PROTECT(N = AS_INTEGER(N));
  PROTECT(M = AS_INTEGER(M));
  PROTECT(alpha = AS_NUMERIC(alpha));
  PROTECT(mu = AS_NUMERIC(mu));
  n = *INTEGER(N);
  m = *INTEGER(M); 
    if(m>1){
   PROTECT(X = allocMatrix(REALSXP, n+1, m));
   PROTECT(Y = allocMatrix(REALSXP, n+1, m));
   PROTECT(Z = allocMatrix(REALSXP, n+1, m));
   PROTECT(Val = allocMatrix(REALSXP, n+1, 3*m));
        }
  else{
   PROTECT(X = NEW_NUMERIC(n+1)); 
   PROTECT(Y = NEW_NUMERIC(n+1)); 
   PROTECT(Z = NEW_NUMERIC(n+1));
   PROTECT(Val = allocMatrix(REALSXP, n+1, 3));
   }
  rX = REAL(X);
  rY = REAL(Y);
  rZ = REAL(Z);
  rVal = REAL(Val);
  rx0 = REAL(x0); 
  ry0 = REAL(y0); 
  rz0 = REAL(z0);
  Alpha = *REAL(alpha);
  Mu = *REAL(mu);  
  for(j=0; j<m; j++){
   rX[j*(n+1)] = rx0[j]; 
   rY[j*(n+1)] = ry0[j];
   rZ[j*(n+1)] = rz0[j];
   rVal[(n+1)*j] = rx0[j];
   rVal[(n+1)*(j+m)] = ry0[j];
   rVal[(n+1)*(j+2*m)] = rz0[j];
   }     
  T = *REAL(t0);
  DELTA = *REAL(delta);
  sd = sqrt(DELTA);
  GetRNGstate();
    for(i=1; i< n+1; i++){
	 T = T + DELTA; 
     for(j=0; j< m; j++){
	  Z1 = rnorm(0,sd);
	  Z2 = rnorm(0,sd); 
	  Z3 = rnorm(0,sd);
	  itrx = rX[i + j*(n+1) - 1]; 
	  itry = rY[i + j*(n+1) - 1]; 
	  itrz = rZ[i + j*(n+1) - 1]; 
	  a1 = fevalxyz(T,itrx,itry,itrz,A1, rho);
	  s1 = fevalxyz(T,itrx,itry,itrz,S1, rho);
	  s11 = fevalxyz(T+DELTA,itrx+a1*DELTA+s1*Z1,itry,itrz,S1,rho);
	  ss11 = fevalxyz(T+DELTA,itrx+a1*DELTA+s1*Z1,itry,itrz,A1,rho)- Mu *fevalxyz(T+DELTA,itrx+a1*DELTA+s1*Z1,itry,itrz,S1,rho) * fevalxyz(T+DELTA,itrx+a1*DELTA+s1*Z1,itry,itrz,S1x,rho);
	  ss12 = fevalxyz(T,itrx,itry,itrz,A1,rho) - Mu *fevalxyz(T,itrx,itry,itrz,S1,rho) * fevalxyz(T,itrx,itry,itrz,S1x,rho);	  
	  a2 = fevalxyz(T,itrx,itry,itrz,A2,rho);
	  s2 = fevalxyz(T,itrx,itry,itrz,S2,rho);
	  s21 = fevalxyz(T+DELTA,itrx,itry+a2*DELTA+s2*Z2,itrz,S2,rho);
	  ss21 = fevalxyz(T+DELTA,itrx,itry+a2*DELTA+s2*Z2,itrz,A2,rho)- Mu *fevalxyz(T+DELTA,itrx,itry+a2*DELTA+s2*Z2,itrz,S2,rho) * fevalxyz(T+DELTA,itrx,itry+a2*DELTA+s2*Z2,itrz,S2x,rho);
	  ss22 = fevalxyz(T,itrx,itry,itrz,A2,rho) - Mu *fevalxyz(T,itrx,itry,itrz,S2,rho) * fevalxyz(T,itrx,itry,itrz,S2x,rho);	  
	  a3 = fevalxyz(T,itrx,itry,itrz,A3,rho);
	  s3 = fevalxyz(T,itrx,itry,itrz,S3,rho);
	  s31 = fevalxyz(T+DELTA,itrx,itry,itrz+a3*DELTA+s3*Z3,S3,rho);
	  ss31 = fevalxyz(T+DELTA,itrx,itry,itrz+a3*DELTA+s3*Z3,A3,rho)- Mu *fevalxyz(T+DELTA,itrx,itry,itrz+a3*DELTA+s3*Z3,S3,rho) * fevalxyz(T+DELTA,itrx,itry,itrz+a3*DELTA+s3*Z3,S3x,rho);
	  ss32 = fevalxyz(T,itrx,itry,itrz,A3,rho) - Mu *fevalxyz(T,itrx,itry,itrz,S3,rho) * fevalxyz(T,itrx,itry,itrz,S3x,rho);	  
      rX[i + (n+1)*j] = itrx + (Alpha*ss11+(1-Alpha)*ss12)*DELTA + (Mu*s11+(1-Mu)*s1)*Z1;
      rY[i + (n+1)*j] = itry + (Alpha*ss21+(1-Alpha)*ss22)*DELTA + (Mu*s21+(1-Mu)*s2)*Z2;
      rZ[i + (n+1)*j] = itrz + (Alpha*ss31+(1-Alpha)*ss32)*DELTA + (Mu*s31+(1-Mu)*s3)*Z3;  
	  rVal[i+(n+1)*j] = rX[i+(n+1)*j];
	  rVal[i+(n+1)*(j+m)] = rY[i+(n+1)*j];
	  rVal[i+(n+1)*(j+2*m)] = rZ[i+(n+1)*j];	  
	  }
	}
  PutRNGstate();
  UNPROTECT(13);
  return Val;
}

     
/** Accessor functions. For internal use only
 * Evaluating R expressions from C
 * Adapted from the code in man/doc/R-exts (p. 118)  
 * fevalx
 * @param t time variable
 * @paramx  space variable
 * @param f a SEXP to a R function
 * @param rho the environment `f' is going to be evaluated
 * @return value of f(t,x)
 * @author S.M. Iacus in pacakge sde (2009).
 * @Ref S.M. Iacus: Simulation and Inference for Stochastic Differential Equations With R Examples, ISBN 978-0-387-75838-1, Springer, NY. 
*/

double fevalx(double t, double x, SEXP f, SEXP rho)
{
    double val= 0.0;
    SEXP R_fcall, tpar, xpar; 
	
	PROTECT(tpar = allocVector(REALSXP, 1));
	PROTECT(xpar = allocVector(REALSXP, 1));
    REAL(tpar)[0] = t;
    REAL(xpar)[0] = x;
   
	PROTECT(R_fcall = allocList(3));
	SETCAR(R_fcall, f);
	SET_TYPEOF(R_fcall, LANGSXP);

	SETCADR(R_fcall, tpar);
	SETCADDR(R_fcall, xpar);
    val = *REAL(eval(R_fcall, rho));
    UNPROTECT(3);
		
    return(val);

}

/** Accessor functions. For internal use only
 * Evaluating R expressions from C
 * Adapted from the code in man/doc/R-exts (p. 118)  
 * fevalx
 * @param t time variable
 * @param x,y  space variable
 * @param f a SEXP to a R function
 * @param rho the environment `f' is going to be evaluated
 * @return value of f(t,x,y)
 * @author A.C. Guidoum  
*/

double fevalxy(double t, double x, double y, SEXP f, SEXP rho)
{
    double val= 0.0;
    SEXP R_fcall, tpar, xpar, ypar; 
	
	PROTECT(tpar = allocVector(REALSXP, 1));
	PROTECT(xpar = allocVector(REALSXP, 1));
	PROTECT(ypar = allocVector(REALSXP, 1));
    REAL(tpar)[0] = t;
    REAL(xpar)[0] = x;
    REAL(ypar)[0] = y;
   
	PROTECT(R_fcall = allocList(4));
	SETCAR(R_fcall, f);
	SET_TYPEOF(R_fcall, LANGSXP);

	SETCADR(R_fcall, tpar);
	SETCADDR(R_fcall, xpar);
	SETCADDDR(R_fcall, ypar);
    val = *REAL(eval(R_fcall, rho));
    UNPROTECT(4);
		
    return(val);
}

/** Accessor functions. For internal use only
 * Evaluating R expressions from C
 * Adapted from the code in man/doc/R-exts (p. 118)  
 * fevalx
 * @param t time variable
 * @param x,y,z  space variable
 * @param f a SEXP to a R function
 * @param rho the environment `f' is going to be evaluated
 * @return value of f(t,x,y,z)
 * @author A.C. Guidoum  
*/

double fevalxyz(double t, double x, double y, double z, SEXP f, SEXP rho)
{
    double val= 0.0;
    SEXP R_fcall, tpar, xpar, ypar, zpar; 
	
	PROTECT(tpar = allocVector(REALSXP, 1));
	PROTECT(xpar = allocVector(REALSXP, 1));
	PROTECT(ypar = allocVector(REALSXP, 1));
	PROTECT(zpar = allocVector(REALSXP, 1));
    REAL(tpar)[0] = t;
    REAL(xpar)[0] = x;
    REAL(ypar)[0] = y;
    REAL(zpar)[0] = z;
   
	PROTECT(R_fcall = allocList(5));
	SETCAR(R_fcall, f);
	SET_TYPEOF(R_fcall, LANGSXP);

	SETCADR(R_fcall, tpar);
	SETCADDR(R_fcall, xpar);
	SETCADDDR(R_fcall, ypar);
	SETCAD4R(R_fcall, zpar);
    val = *REAL(eval(R_fcall, rho));
    UNPROTECT(5);
		
    return(val);
}
