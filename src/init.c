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
 * @file   init.c
 * @author A.C. Guidoum and K. Boukhetala
 * @date   2011-2014
 * @brief  Init File for All C Files. 
 */ 

#include <R_ext/Rdynload.h>
#include "Sim.DiffProc.h"

/**
 * Registering native routines (in man/doc/R-exts, p. 95)
*/


static R_CMethodDef cMethods[] = {
   {"Euler1d", (DL_FUNC)& Euler1d, 8},
   {"Euler2d", (DL_FUNC)& Euler2d, 11},
   {"Euler3d", (DL_FUNC)& Euler3d, 14},
   {"Milstein1d", (DL_FUNC)& Milstein1d, 9},
   {"Milstein2d", (DL_FUNC)& Milstein2d, 13},
   {"Milstein3d", (DL_FUNC)& Milstein3d, 17},
   {"SMilstein1d", (DL_FUNC)& SMilstein1d, 12},
   {"SMilstein2d", (DL_FUNC)& SMilstein2d, 19},
   {"SMilstein3d", (DL_FUNC)& SMilstein3d, 26},
   {"Sts1d", (DL_FUNC)& Sts1d, 14},
   {"Sts2d", (DL_FUNC)& Sts2d, 23},
   {"Sts3d", (DL_FUNC)& Sts3d, 32},
   {"Heun1d", (DL_FUNC)& Heun1d, 8},
   {"Heun2d", (DL_FUNC)& Heun2d, 11},
   {"Heun3d", (DL_FUNC)& Heun3d, 14},
   {"Rk1d", (DL_FUNC)& Rk1d, 9},
   {"Rk2d", (DL_FUNC)& Rk2d, 12},
   {"Rk2d", (DL_FUNC)& Rk3d, 15},
   {"Predcorr1d", (DL_FUNC)& Predcorr1d, 11},
   {"Predcorr2d", (DL_FUNC)& Predcorr2d, 15},
   {"Predcorr3d", (DL_FUNC)& Predcorr3d, 19},
   {NULL, NULL, 0},
};


// void R_init_ifs(DllInfo *info)
// {
    // R_registerRoutines(info, R_CDef, NULL, NULL, NULL);
// }

// void R_init_markovchain(DllInfo *info) {
  // R_registerRoutines(info, R_CDef, NULL, NULL, NULL);
  // R_useDynamicSymbols(info, FALSE);
  // R_forceSymbols(info, TRUE); 
// }

void R_init_markovchain(DllInfo *info) {
  R_registerRoutines(info, cMethods, NULL, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
  R_forceSymbols(info, TRUE);
}

