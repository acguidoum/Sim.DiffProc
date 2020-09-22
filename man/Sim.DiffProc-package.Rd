
\newcommand{\HR}{\href{http://www.usthb.dz/en/}}
\newcommand{\RR}{\href{http://www.r-project.org}}

\name{Sim.DiffProc-package}
\alias{Sim.DiffProc-package}
\alias{Sim.DiffProc}
\docType{package}
\title{
Simulation of Diffusion Processes
}
\description{
A package for symbolic and numerical computations on scalar and multivariate systems of stochastic differential equations. It provides users with a wide range of tools to simulate, estimate, analyze, and visualize the dynamics of these systems in both forms Itô and Stratonovich. Statistical analysis with parallel Monte-Carlo and moment equations methods of SDE's. Enabled many searchers in different domains to use these equations to modeling practical problems in financial and actuarial modeling and other areas of application, e.g., modeling and simulate of first passage time problem in shallow water using the attractive center (Boukhetala K, 1996) ISBN:1-56252-342-2. 
}
\details{
\tabular{ll}{
Package: \tab Sim.DiffProc \cr
Type: \tab Package \cr
Version: \tab 4.7 \cr
Date: \tab 2020-09-22 \cr
License: \tab GPL (>= 2) \cr
Depends: \tab R (>= 3.0.0) \cr
Imports: \tab Deriv (>= 3.8.0), MASS (>= 7.3-30), parallel \cr
Suggests: \tab deSolve (>= 1.11), knitr (>= 1.10), rgl (>= 0.93.991), rmarkdown (>= 0.8), scatterplot3d (>= 0.3-36), sm (>= 2.2-5.3) \cr
Classification/MSC:\tab 37H10, 37M10, 60H05, 60H10, 60H35, 60J60, 65C05, 68N15, 68Q10 \cr
}

There are main types of functions in this package:
  \enumerate{
    \item Simulation of solution to 1,2 and 3-dim stochastic differential equations of Itô and Stratonovich types, with different methods.
    \item Simulation of solution to 1,2 and 3-dim diffusion bridge of Itô and Stratonovich types, with different methods. 
	\item Simulation the first-passage-time (f.p.t) in 1,2 and 3-dim sde of Itô and Stratonovich types.
	\item Calculate symbolic ODE's of moment equations (means and variances-covariance) for 1,2 and 3-dim SDE's.
	\item Monte-Carlo replicates of a statistic applied to 1,2 and 3-dim SDE's at any time t.
    \item Computing the basic statistics (mean, var, median, ...) of the processes at any time t using the Monte Carlo method.
    \item Random number generators (RN's) to generate 1,2 and 3-dim sde of Itô and Stratonovich types.
	\item Approximate the transition density  1,2 and 3-dim of the processes at any time t.
	\item Approximate the density of first-passage-time in 1,2 and 3-dim SDE's.
    \item Computing the Itô and Stratonovich stochastic integrals.
    \item Estimate drift and diffusion parameters by the method of maximum pseudo-likelihood of the one-dim stochastic differential equation.	
	\item Converting Sim.DiffProc objects to LaTeX.
    \item Displaying an object inheriting from class \code{"sde"} (1,2 and 3 dim).
  }

For more examples see \code{demo(Sim.DiffProc)}, and for an overview of this package, see \code{\link[utils]{browseVignettes}('Sim.DiffProc')} for more informations. 
}
\section{Requirements}{
\R version >= 3.0.0
}
\author{
A.C. Guidoum \email{acguidoum@usthb.dz} and K. Boukhetala \email{kboukhetala@usthb.dz} (Dept. Probability and Statistics, USTHB, Algeria).\cr 
Please send comments, error reports, etc. to the author via the addresses email.
}
\section{Licence}{
This package and its documentation are usable under the terms of the "GNU General Public License", a copy of which is distributed with the package.
}
\references{
   Argyrakisa, P. and G.H. Weiss (2006).
   A first-passage time problem for many random walkers.
   \emph{Physica A}. \bold{363}, 343--347.
  
   Aytug H., G. J. Koehler (2000). 
   New stopping criterion for genetic algorithms. 
   \emph{European Journal of Operational Research}, \bold{126}, 662--674.
   
  
   Boukhetala, K. (1994). 
   Simulation study of a dispersion about an attractive centre. 
   In proceedings of 11th Symposium Computational Statistics, edited by R.Dutter and W.Grossman, Wien , Austria, 128--130.  
  
   Boukhetala, K. (1996).
   Modelling and simulation of a dispersion pollutant with attractive centre.
   ed by Computational Mechanics Publications, Southampton ,U.K and Computational Mechanics Inc, Boston, USA, 245--252.
  
   Boukhetala, K. (1998a). 
   Estimation of the first passage time distribution for a simulated diffusion process.
   \emph{Maghreb Math.Rev}, \bold{7}(1), 1--25.
  
   Boukhetala, K. (1998b). 
   Kernel density of the exit time in a simulated diffusion. 
   \emph{les Annales Maghrebines De L ingenieur}, \bold{12}, 587--589.
   
   Ding, M. and G. Rangarajan. (2004). 
   First Passage Time Problem: A Fokker-Planck Approach. 
   \emph{New Directions in Statistical Physics}. ed by L. T. Wille. Springer. 31--46.
   
   Ait-Sahalia, Y. (1999). 
   Transition densities for interest rate and other nonlinear diffusions.
   \emph{The Journal of Finance}, 54, 1361--1395.
   
   Ait-Sahalia, Y. (2002). 
   Maximum likelihood estimation of discretely sampled diffusions: a closed-form approximation approach. 
   \emph{Econometrica}. 70, 223--262.
  
   Roman, R.P., Serrano, J. J., Torres, F. (2008).
   First-passage-time location function: Application to determine first-passage-time densities in diffusion processes.
   \emph{Computational Statistics and Data Analysis}. \bold{52}, 4132--4146.
  
   Roman, R.P., Serrano, J. J., Torres, F. (2012). 
   An R package for an efficient approximation of first-passage-time densities for diffusion processes based on the FPTL function. 
   \emph{Applied Mathematics and Computation}, \bold{218}, 8408--8428.
    
   Kessler, M. (1997).  
   Estimation of an ergodic diffusion from discrete observations. 
   \emph{Scand. J. Statist.}, 24, 211-229.
  
   Gardiner, C. W. (1997). 
   \emph{Handbook of Stochastic Methods}.
   Springer-Verlag, New York.
  
   Friedman, A. (1975). 
   \emph{Stochastic differential equations and applications}.
   Volume 1, ACADEMIC PRESS.
  
   Henderson, D. and Plaschko,P. (2006). 
   \emph{Stochastic differential equations in science and engineering}.
   World Scientific.
   
   Croissant, Y. (2014). 
   Ecdat: Data sets for econometrics. 
   \emph{\R package version 0.2-5}.
   
   Vasicek, O. (1977). 
   An Equilibrium Characterization of the Term Structure.
   \emph{Journal of Financial Economics},  5, 177--188. 
  
   Allen, E. (2007).
   \emph{Modeling with Itô stochastic differential equations}. 
   Springer-Verlag.

   Jedrzejewski, F. (2009). 
   \emph{Modeles aleatoires et physique probabiliste}. 
   Springer-Verlag.
  
   Iacus, S.M. (2008). 
   \emph{Simulation and inference for stochastic differential equations: with R examples}. 
   Springer-Verlag, New York.
   
   Iacus, S.M. (2014). 
   sde: Simulation and Inference for Stochastic Differential Equations. 
   \emph{\R package version 2.0.13}. 
   
   Brouste, A. et al. (2014).
   The yuima Project: A Computational Framework for Simulation and Inference of Stochastic Differential Equations.
  \emph{Journal of Statistical Software}, \bold{57}(4).
  
   Kloeden, P.E, and Platen, E. (1989). 
   A survey of numerical methods for stochastic differential equations. 
   \emph{Stochastic Hydrology and Hydraulics}, \bold{3}, 155--178.
      

   Kloeden, P.E, and Platen, E. (1991a).
   Relations between multiple Itô and stratonovich integrals.
   \emph{Stochastic Analysis and Applications}, \bold{9}(3), 311--321.

   Kloeden, P.E, and Platen, E. (1991b).
   Stratonovich and Itô stochastic taylor expansions. 
   \emph{Mathematische Nachrichten}, \bold{151}, 33--50.

   Kloeden, P.E, and Platen, E. (1995).
   \emph{Numerical Solution of Stochastic Differential Equations}. 
   Springer-Verlag, New York.

   Oksendal, B. (2000).
   \emph{Stochastic Differential Equations: An Introduction with Applications}.
   5th edn. Springer-Verlag, Berlin.
   
   B.L.S. Prakasa Rao. (1999).
   \emph{Statistical Inference for Diffusion Type Processes}.
   Arnold, London and Oxford University press, New York.

   Kutoyants, Y.A. (2004).
   \emph{Statistical Inference for Ergodic Diffusion Processes}.
   Springer, London.
   
   Sorensen, H. (2000).
   Inference for Diffusion Processes and Stochastic Volatility Models.
   Ph.D. thesis, Department of Mathematical Sciences, University of Copenhagen.
   
   Sorensen, H. (2002).
   Estimation of diffusion parameters for discretely observed diffusion processes.
   \emph{Bernoulli}, 8, 491--508.

   Sorensen, H. (2004).
   Parametric inference for diffusion processes observed at discrete points in time: a survey.
   \emph{International Statistical Review}, 72, 337--354.

   Platen, E. (1980).
   Weak convergence of approximations of Itô integral equations. 
   \emph{Z Angew Math Mech}. \bold{60}, 609--614.

   Platen, E. and Bruti-Liberati, N. (2010).
   \emph{Numerical Solution of Stochastic Differential Equations with Jumps in Finance.} 
   Springer-Verlag, New York.

   Itô, Y, and Mitsui, T. (1993).
   Simulation of Stochastic Differential Equations. 
   \emph{The Annals of the Institute of Statistical Mathematics}, \bold{3}, 419--432.
   
   Risken, H. (2001).
   \emph{The Fokker Planck Equation : Methods of Solutions and Applications}. 
   2nd edition, Springer Series in Synergetics.
   
  Dacunha, D.C. and Florens, D.Z. (1986). 
  Estimation of the Coefficients of a Diffusion from Discrete Observations. 
  \emph{Stochastics}. 19, 263--284.  

  Dohnal, G. (1987).
  On estimating the diffusion coefficient.
  \emph{J. Appl.Prob.}, 24, 105--114.

  Genon, V.C. (1990).
  Maximum constrast estimation for diffusion processes from discrete observation.
  \emph{Statistics}, 21, 99--116.
   
   Protter, P. (2005). 
   \emph{Stochastic Integration and Differential Equations.}
   2nd edn. Springer-Verlag, New York.
   
   Bladt, M. and Sorensen, M. (2007).
   Simple simulation of diffusion bridges with application to likelihood inference for diffusions.
   \emph{Working Paper, University of Copenhagen}. 
   
   Ozaki, T. (1992). 
   A bridge between nonlinear time series models and nonlinear stochastic dynamical systems: A local linearization approach. 
   \emph{Statistica Sinica}, 2, 25-83.
   
   Shoji, L., Ozaki, T. (1998). 
   Estimation for nonlinear stochastic differential equations by a local linearization method. 
   \emph{Stochastic Analysis and Applications}, 16, 733-752.
   
   Nicolau, J. (2004).
   Introduction to the estimation of stochastic differential equations based on discrete observations.
   \emph{Autumn School and International Conference, Stochastic Finance}.
   
   F C Klebaner, F.C. (2005).
   \emph{Introduction to stochastic calculus with application.} 
   2nd edn. Imperial College Press (ICP).
     
  Henderson, D. and Plaschko, P. (2006). 
  \emph{Stochastic differential equations in science and engineering.}
  World Scientific.
}
\keyword{package}
\seealso{
sde, yumia, DiffusionRgqd, DiffusionRjgqd, DiffusionRimp, QPot, diffeqr, fptdApprox.
}

