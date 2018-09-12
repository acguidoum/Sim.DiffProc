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
 * @file   Sim.DiffProc.h
 * @author A.C. Guidoum and K. Boukhetala
 * @date   2011-2014
 * @brief  Header File for Numerical Methods for Multidimensional Stochastic Differential Equations. 
 */
 
/** 
 * @Functions accessed from .Call() 
 */

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <math.h>
#include "locale.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

# define max(a,b) ((a) > (b) ? (a) : (b))
# define min(a,b) ((a) < (b) ? (a) : (b))


# define abs9(a) (a > 0 ? a:-a)


SEXP Euler1d(SEXP x0, SEXP t0, SEXP delta, SEXP N, SEXP M,
             SEXP A, SEXP S, SEXP rho);

SEXP Euler2d(SEXP x0, SEXP y0, SEXP t0, SEXP delta, SEXP N, SEXP M,
             SEXP A1, SEXP S1, SEXP A2, SEXP S2, SEXP rho);	

SEXP Euler3d(SEXP x0, SEXP y0, SEXP z0, SEXP t0, SEXP delta, SEXP N, 
             SEXP M, SEXP A1, SEXP S1, SEXP A2, SEXP S2, SEXP A3, 
			 SEXP S3,SEXP rho);			 
			 				   
SEXP Milstein1d(SEXP x0, SEXP t0, SEXP delta, SEXP N, SEXP M,
                SEXP A, SEXP S, SEXP Sx, SEXP rho);	
				
SEXP Milstein2d(SEXP x0, SEXP y0, SEXP t0, SEXP delta, SEXP N, 
                SEXP M, SEXP A1, SEXP S1, SEXP S1x, SEXP A2, 
				SEXP S2, SEXP S2x, SEXP rho);
				
SEXP Milstein3d(SEXP x0, SEXP y0, SEXP z0, SEXP t0, SEXP delta, SEXP N, 
                SEXP M, SEXP A1, SEXP S1, SEXP S1x, SEXP A2, 
				SEXP S2, SEXP S2x, SEXP A3, SEXP S3, SEXP S3x, SEXP rho);

SEXP SMilstein1d(SEXP x0, SEXP t0, SEXP delta, SEXP N, SEXP M,
                 SEXP A, SEXP Ax, SEXP Axx, SEXP S, SEXP Sx, 
				 SEXP Sxx, SEXP rho);
				 
SEXP SMilstein2d(SEXP x0, SEXP y0, SEXP t0, SEXP delta, SEXP N, 
                 SEXP M, SEXP A1, SEXP A1x, SEXP A1xx, SEXP A2, 
				 SEXP A2x, SEXP A2xx, SEXP S1, SEXP S1x, SEXP S1xx, 
				 SEXP S2, SEXP S2x, SEXP S2xx, SEXP rho);

SEXP SMilstein3d(SEXP x0, SEXP y0, SEXP z0, SEXP t0, SEXP delta, SEXP N, 
                 SEXP M, SEXP A1, SEXP A1x, SEXP A1xx, SEXP A2, 
				 SEXP A2x, SEXP A2xx, SEXP A3, SEXP A3x, SEXP A3xx, SEXP S1, 
				 SEXP S1x, SEXP S1xx, SEXP S2, SEXP S2x, SEXP S2xx, 
				 SEXP S3, SEXP S3x, SEXP S3xx, SEXP rho);				 

SEXP Sts1d(SEXP x0, SEXP t0, SEXP delta, SEXP N, SEXP M,
           SEXP A, SEXP Ax, SEXP Axx, SEXP S, SEXP Sx, 
		   SEXP Sxx, SEXP Z, SEXP U, SEXP rho);
		   
SEXP Sts2d(SEXP x0, SEXP y0, SEXP t0, SEXP delta, SEXP N, SEXP M, 
           SEXP A1, SEXP A1x, SEXP A1xx, SEXP A2, SEXP A2x, SEXP A2xx, 
		   SEXP S1, SEXP S1x, SEXP S1xx, SEXP S2, SEXP S2x, SEXP S2xx, 
		   SEXP Z1, SEXP U1, SEXP Z2, SEXP U2, SEXP rho);	

SEXP Sts3d(SEXP x0, SEXP y0, SEXP z0, SEXP t0, SEXP delta, SEXP N, 
           SEXP M, SEXP A1, SEXP A1x, SEXP A1xx, SEXP A2, 
		   SEXP A2x, SEXP A2xx, SEXP A3, SEXP A3x, SEXP A3xx, SEXP S1, 
		   SEXP S1x, SEXP S1xx, SEXP S2, SEXP S2x, SEXP S2xx, 
		   SEXP S3, SEXP S3x, SEXP S3xx, SEXP Z1, SEXP U1, SEXP Z2, 
		   SEXP U2, SEXP Z3, SEXP U3, SEXP rho);		   

SEXP Heun1d(SEXP x0, SEXP t0, SEXP delta, SEXP N, SEXP M,
            SEXP A, SEXP S, SEXP rho);	
			
SEXP Heun2d(SEXP x0, SEXP y0, SEXP t0, SEXP delta, SEXP N, SEXP M,
            SEXP A1, SEXP S1, SEXP A2, SEXP S2, SEXP rho);	

SEXP Heun3d(SEXP x0, SEXP y0, SEXP z0, SEXP t0, SEXP delta, SEXP N, 
            SEXP M, SEXP A1, SEXP S1, SEXP A2, SEXP S2, SEXP A3, 
			SEXP S3, SEXP rho);			

SEXP Rk1d(SEXP x0, SEXP t0, SEXP delta, SEXP N, SEXP M,
          SEXP A, SEXP S, SEXP Order, SEXP rho);

SEXP Rk2d(SEXP x0, SEXP y0, SEXP t0, SEXP delta, SEXP N, SEXP M,
          SEXP A1, SEXP S1, SEXP A2, SEXP S2, SEXP Order, SEXP rho);

SEXP Rk3d(SEXP x0, SEXP y0, SEXP z0, SEXP t0, SEXP delta, SEXP N, 
          SEXP M, SEXP A1, SEXP S1, SEXP A2, SEXP S2, SEXP A3, 
	      SEXP S3, SEXP Order, SEXP rho);		  
			 
SEXP Predcorr1d(SEXP x0, SEXP t0, SEXP delta, SEXP N, SEXP M,
                SEXP alpha, SEXP mu, SEXP A, SEXP S, SEXP SS, SEXP rho);	

SEXP Predcorr2d(SEXP x0, SEXP y0, SEXP t0, SEXP delta, SEXP N, SEXP M,
                SEXP alpha, SEXP mu, SEXP A1, SEXP S1, SEXP S1x,  
				SEXP A2, SEXP S2, SEXP S2x, SEXP rho);

SEXP Predcorr3d(SEXP x0, SEXP y0, SEXP z0, SEXP t0, SEXP delta, SEXP N, 
                SEXP M, SEXP alpha, SEXP mu, SEXP A1, SEXP S1, SEXP S1x, 
				SEXP A2, SEXP S2, SEXP S2x, SEXP A3, SEXP S3, SEXP S3x, 
				SEXP rho);				

double fevalx(double t, double x, SEXP f, SEXP rho);

double fevalxy(double t, double x, double y, SEXP f, SEXP rho);

double fevalxyz(double t, double x, double y, double z, SEXP f, SEXP rho);
 
 
 
