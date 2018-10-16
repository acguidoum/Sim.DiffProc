## Sun Oct 14 00:31:28 2018
## Original file Copyright © 2018 A.C. Guidoum, K. Boukhetala
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


#############
#############
## .noGenerics <- TRUE

.onLoad <- function(libname, pkgname)
           {
    if (!interactive()) options(rgl.useNULL = TRUE)
 }

## .onUnload <- function(libpath) {
##     library.dynam.unload("Sim.DiffProc", libpath)
## }

# .onAttach <- function(libname, pkgname) {
    # packageStartupMessage(paste0("This is package 'Sim.DiffProc', v",packageVersion(pkgname)));
# }

.onAttach <- function(library, pkg) {
    packageStartupMessage("Package 'Sim.DiffProc', version 4.2\nbrowseVignettes('Sim.DiffProc') for more informations.")
	invisible()
}
