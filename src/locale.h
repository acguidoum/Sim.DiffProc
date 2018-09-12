/** 
 * Sat Aug 09 02:27:46 2014
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
 * @file   locale.h
 * @author A.C. Guidoum and K. Boukhetala
 * @date   2011-2014
 * @brief  header file for error messages. 
 */
 



/** 
 * @Localization 
 */
#include <R.h>
#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("Sim.DiffProc", String)
#else
#define _(String) (String)
#endif

