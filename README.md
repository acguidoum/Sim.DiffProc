<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Project Status](https://www.repostatus.org/badges/latest/active.svg?style=popout)](https://cran.r-project.org/src/contrib/Archive/Sim.DiffProc/)
[![Travis Build Status](https://travis-ci.org/acguidoum/Sim.DiffProc.svg?branch=master)](https://travis-ci.org/acguidoum/Sim.DiffProc)
[![Build status](https://ci.appveyor.com/api/projects/status/16a70vyf8rk7nn1i?svg=true)](https://ci.appveyor.com/project/acguidoum/sim-diffproc-xal8n) [![codecov](https://codecov.io/gh/acguidoum/Sim.DiffProc/branch/master/graph/badge.svg)](https://codecov.io/gh/acguidoum/Sim.DiffProc)


------------------------------------------------------------------------

[![CRAN](https://img.shields.io/cran/l/devtools.svg?style=popout)](https://cran.r-project.org/web/licenses/GPL-2)
[![minimal R version](https://img.shields.io/badge/R%3E%3D-2.15.1-blue.svg?style=flat-plastic)](https://cran.r-project.org/) 

------------------------------------------------------------------------
![Github](https://img.shields.io/badge/Github-4.2 (2018.10.15)-blue.svg)
![](http://www.r-pkg.org/badges/version-last-release/Sim.DiffProc?color=blue) [![Rdoc](http://www.rdocumentation.org/badges/version/Sim.DiffProc)](http://www.rdocumentation.org/packages/Sim.DiffProc)


![](https://cranlogs.r-pkg.org/badges/grand-total/Sim.DiffProc?color=yellow)
![](https://cranlogs.r-pkg.org/badges/Sim.DiffProc?color=yellow)
![](https://cranlogs.r-pkg.org/badges/last-week/Sim.DiffProc?color=yellow)
![](https://cranlogs.r-pkg.org/badges/last-day/Sim.DiffProc?color=yellow)

------------------------------------------------------------------------


Package Overview
---------------------

The package [Sim.DiffProc](https://cran.r-project.org/package=Sim.DiffProc) is an object created in R for symbolic and numerical computations on scalar and multivariate systems of stochastic differential equations. It provides users with a wide range of tools to simulate, estimate, analyze, and visualize the dynamics of these systems in both forms Ito and Stratonovich. The project was officially launched in September 2010 and is under active development by the authors. The current feature set of the package can be split in more main categories: Computing the stochastic integrals of Ito or Stratonovich type. Simulation sde's and bridge sde's of Ito or Stratonovich type (1,2 and 3-dim), with different methods. Approximate transition density and random number generators for SDE's. Density approximation for First-passage-time (f.p.t) in SDE's (1,2 and 3-dim). Statistical analysis with Parallel Monte-Carlo and moment equations methods of SDE's (1,2 and 3-dim). Estimate drift and diffusion parameters using pseudo-maximum likelihood estimators of 1-dim SDE's. Displaying an object inheriting from a class of SDE's.

The package includes the following categories (where `k=1,2,3`):

1. [`snssdekd()` & `dsdekd()` & `rsdekd()`- Monte-Carlo Simulation and Analysis of Stochastic Differential Equations](https://CRAN.R-project.org/package=Sim.DiffProc/vignettes/snssde.html).
2. [`bridgesdekd()` & `dsdekd()` & `rsdekd()` - Constructs and Analysis of Bridges Stochastic Differential Equations](https://CRAN.R-project.org/package=Sim.DiffProc/vignettes/bridgesde.html).
3. [`fptsdekd()` & `dfptsdekd()` - Monte-Carlo Simulation and Kernel Density Estimation of First passage time](https://CRAN.R-project.org/package=Sim.DiffProc/vignettes/fptsde.html).
4. [`MCM.sde()` & `MEM.sde()` - Parallel Monte-Carlo and Moment Equations for SDEs](https://CRAN.R-project.org/package=Sim.DiffProc/vignettes/mcmsde.html).
5. [`TEX.sde()` - Converting Sim.DiffProc Objects to LaTeX](https://CRAN.R-project.org/package=Sim.DiffProc/vignettes/sdetotex.html).
6. [`fitsde()` - Parametric Estimation of 1-D Stochastic Differential Equation](https://CRAN.R-project.org/package=Sim.DiffProc/vignettes/fitsde.html).

Obtaining and installation
-----------------------

As `Sim.DiffProc` is an `R` package, it requires `R version 2.15.1` or higher to be installed, distributed as open source software under the GPL-2/GPL-3 license. The package is available from CRAN at URL https://CRAN.R-project.org/package=Sim.DiffProc, or from GitHub at URL https://github.com/acguidoum/Sim.DiffProc. To download, install and load the current release, just type the code below in your current `R` session:

```r
install.packages("Sim.DiffProc")
## Or 
install.packages("devtools")
devtools::install_github("acguidoum/Sim.DiffProc")
library("Sim.DiffProc")
```

Documentation and Examples
--------------------------

It is a requirement of the R packaging system that every function and data set in a package has a help page. The [Sim.DiffProc](https://cran.r-project.org/package=Sim.DiffProc) package follows this  requirement strictly. In addition to the help pages, the package includes vignettes and demonstration scripts. First read the package vignette Then read the [reference manual.](https://CRAN.R-project.org/package=Sim.DiffProc/Sim.DiffProc.pdf)



```r
browseVignettes(package = "Sim.DiffProc")
```
and 

```r
demo(package = "Sim.DiffProc")
```

Collaboration and citation
-----

Obviously, the package leaves many other fields of stochastic modeling with Ito and Stratonovich SDE's untouched. For this situation to change, we hope that experts in their field will join their efforts to ours and contribute code to the Sim.DiffProc project. The project will continue to grow and improve by the authors to the community of developers and users. If you use [Sim.DiffProc](https://cran.r-project.org/package=Sim.DiffProc) please cite the software in publications;
use `citation()` for information on how to cite the software;

```r
citation("Sim.DiffProc")
 
To cite package 'Sim.DiffProc' in publications use:

  Arsalane Chouaib Guidoum and Kamal Boukhetala (2018).
  Sim.DiffProc: Simulation of Diffusion Processes.R package
  version 4.1.https://cran.r-project.org/package=Sim.DiffProc

A BibTeX entry for LaTeX users is

  @Manual{,
    title = {Sim.DiffProc: Simulation of Diffusion Processes.},
    author = {Arsalane Chouaib Guidoum and Kamal Boukhetala},
    year = {2018},
    note = {R package version 4.2},
    url = {https://cran.r-project.org/package=Sim.DiffProc},
  }
```

Note
----

Please send comments, error reports, etc. to the author via the addresses email.

References
----------

1. Guidoum AC, Boukhetala K (2017). Performing Parallel Monte Carlo and Moment Equations Methods for Ito and Stratonovich Stochastic Differential Systems: R Package Sim.DiffProc. Preprint submitted to Journal of Statistical Software.

2. Guidoum AC, Boukhetala K (2018). Sim.DiffProc: Simulation of Diffusion Processes. R package version 4.2, URL https://cran.r-project.org/package=Sim.DiffProc.


