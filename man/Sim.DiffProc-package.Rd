
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
Version: \tab 4.1 \cr
Date: \tab 2018-09-09 \cr
License: \tab GPL (>= 2) \cr
Depends: \tab R (>= 2.15.1) \cr
Imports: \tab Deriv (>= 3.8.0), MASS, parallel, rgl (>= 0.66), scatterplot3d (>= 0.3-36) \cr
Suggests: \tab deSolve (>= 1.11), knitr, sm (>= 2.2-5.3) \cr
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
    \item Computing the stochastic integrals of Itô and Stratonovich types.
    \item Estimate drift and diffusion parameters by the method of maximum pseudo-likelihood of the 1-dim stochastic differential equation.	
	\item Converting Sim.DiffProc objects to LaTeX.
    \item Displaying an object inheriting from class \code{"sde"} (1,2 and 3 dim).
  }
}
\section{Main Features}{
\bold{Stochastic integrals:}\cr\cr
We consider a simple example to simulation Itô integral, used \code{\link{st.int}} function:
\deqn{\int_{t_{0}}^{t} W_{s}^{n} dW_{s} = \frac{1}{n+1} \left[W_{t}^{n+1} - W_{t_{0}}^{n+1}\right]- \frac{n}{2} \int_{t_{0}}^{t} W_{s}^{n-1} ds}{ int(w^n * dw(s),0,T) = w^(n+1) /n+1 - n/2 int(w^(n-1) ds,0,T)}
And the Stratonovich integral
\deqn{\int_{t_{0}}^{t} W_{s}^{n} \circ dW_{s} = \frac{1}{n+1} \left[W_{t}^{n+1} - W_{t_{0}}^{n+1}\right]}{ int(w^n o dw(s),0,T) = w^(n+1) /n+1 }
  \preformatted{
  R> f <- expression( w )
  R> Itô <- st.int(f,type="Ito",M=500,lower=0,upper=1)
  R> Itô
   Itô integral:
        | X(t)   = integral (f(s,w) * dw(s))
        | f(t,w) = w
   Summary:
        | Number of subinterval | N = 1001.
        | Number of simulation  | M = 500.
        | Limits of integration | t in [0,1].
  R> summary(Itô)
   Monte-Carlo Statistics for integral(f(s,w) * dw(s)) at time t = 1
   | f(t,w) = w
                         
  Mean              0.01330
  Variance          0.51102
  Median           -0.28645
  Mode             -0.42772
  First quartile   -0.44666
  Third quartile    0.22534
  Minimum          -0.55198
  Maximum           4.38802
  Skewness          2.27133
  Kurtosis          9.27393
  Coef-variation   53.75783
  3th-order moment  0.82972
  4th-order moment  2.42178
  5th-order moment  7.60355
  6th-order moment 26.72897 
  
  R> str <- st.int(f,type="str",M=500,lower=0,upper=1)
  R> str
  Stratonovich integral:
        | X(t)   = integral (f(s,w) o dw(s))
        | f(t,w) = w
  Summary:
        | Number of subinterval | N = 1001.
        | Number of simulation  | M = 500.
        | Limits of integration | t in [0,1].
  R> summary(str)
   Monte-Carlo Statistics for integral (f(s,w) o dw(s)) at time t = 1
   | f(t,w) = w
                          
   Mean               0.55655
   Variance           0.66663
   Median             0.21223
   Mode               0.08249
   First quartile     0.04269
   Third quartile     0.79322
   Minimum            0.00000
   Maximum            6.70508
   Skewness           2.72896
   Kurtosis          13.34205
   Coef-variation     1.46702
   3th-order moment   1.48532
   4th-order moment   5.92908
   5th-order moment  26.87087
   6th-order moment 138.27901  
  }
\bold{SDE's 1,2 and 3-dim:}\cr\cr
There are thus two widely used types of stochastic calculus, Stratonovich and Itô, differing in respect of the stochastic integral used. Modelling issues typically dictate which version in appropriate, but once one has been chosen a corresponding equation of the other type with the same solutions can be determined. Thus it is possible to switch between the two stochastic calculi. Specifically, the processes \eqn{\{ X_{t}, t \geq 0 \}}{X(t)} solution to the Itô SDE:
\deqn{dX_{t} = f(t,X_{t}) dt + g(t,X_{t}) dW_{t}}{dX(t) = f(t,X(t)) dt + g(t,X(t)) dW(t)}
where \eqn{\{ W_{t}, t \geq 0 \}}{W(t)} is the standard Wiener process or standard Brownian motion, the drift \eqn{f(t,X_{t})}{f(t,X(t))} and diffusion \eqn{g(t,X_{t})}{g(t,X(t))} are known functions that are assumed to be sufficiently regular (Lipschitz, bounded growth) for existence and uniqueness of solution; has the same solutions as the Stratonovich SDE:
\deqn{dX_{t} = \underline{f}(t,X_{t}) dt + g(t,X_{t}) \circ dW_{t}}{dX(t) = bar_f(t,X(t)) dt + g(t,X(t)) o dW(t)}
with the modified drift coefficient
\deqn{\underline{f}(t,X_{t})= f(t,X_{t}) - \frac{1}{2} g(t,X_{t})\frac{\partial g}{\partial x}(t,X_{t})}{bar_f(t,X(t)) = f(t,X(t)) - 1/2 g(t,X(t)) gx(t,X(t))}  
The following examples for different methods of simulation of SDEs (1,2 and 3-dim) use the \code{\link{snssde1d}}, \code{\link{snssde2d}} and \code{\link{snssde3d}} functions.
  \preformatted{
  R> ## 1-dim sde
  R> f <- expression(2*(3-x) )
  R> g <- expression(2*x)
  R> res1 <- snssde1d(drift=f,diffusion=g,M=1000,x0=1)
  R> res1
  Itô Sde 1D:
        | dX(t) = 2 * (3 - X(t)) * dt + 2 * X(t) * dW(t)
  Method:
        | Euler scheme with order 0.5
  Summary:
        | Size of process       | N  = 1001.
        | Number of simulation  | M  = 1000.
        | Initial value         | x0 = 1.
        | Time of process       | t in [0,1].
        | Discretization        | Dt = 0.001.	
  R> res2 <- snssde1d(drift=f,diffusion=g,M=1000,x0=1,type="str")
  R> res2
  Stratonovich Sde 1D:
        | dX(t) = 2 * (3 - X(t)) * dt + 2 * X(t) o dW(t)
  Method:
        | Euler scheme with order 0.5
  Summary:
        | Size of process       | N  = 1001.
        | Number of simulation  | M  = 1000.
        | Initial value         | x0 = 1.
        | Time of process       | t in [0,1].
        | Discretization        | Dt = 0.001.
		
  R> ## 2-dim sde
  R> fx <- expression(x-y, y-x)
  R> gx <- expression(2*y, 2*x)
  R> res2d <- snssde2d(drift=fx,diffusion=gx,x0=c(1,1))
  R> res2d
  Itô Sde 2D:
        | dX(t) = X(t) - Y(t) * dt + 2 * Y(t) * dW1(t)
        | dY(t) = Y(t) - X(t) * dt + 2 * X(t) * dW2(t)
  Method:
        | Euler scheme with order 0.5
  Summary:
        | Size of process       | N  = 1001.
        | Number of simulation  | M  = 1.
        | Initial values        | (x0,y0) = (1,1).
        | Time of process       | t in [0,1].
        | Discretization        | Dt = 0.001.
  R> plot2d(res2d)
		
  R> ## 3-dim sde
  R> fx <- expression(y, 0, 0)
  R> gx <- expression(z, 1, 1)
  R> res3d <- snssde3d(drift=fx,diffusion=gx,M=1000)
  R> res3d
  Itô Sde 3D:
        | dX(t) = Y(t) * dt + Z(t) * dW1(t)
        | dY(t) = 0 * dt + 1 * dW2(t)
        | dZ(t) = 0 * dt + 1 * dW3(t)
  Method:
        | Euler scheme with order 0.5
  Summary:
        | Size of process       | N  = 1001.
        | Number of simulation  | M  = 1000.
        | Initial values        | (x0,y0,z0) = (0,0,0).
        | Time of process       | t in [0,1].
        | Discretization        | Dt = 0.001.
  plot3D(res3d) 
  }
\bold{Bridge SDE's 1,2 and 3-dim:}\cr\cr
Simulation of bridge SDEs (1,2 and 3-dim) with \code{\link{bridgesde1d}}, \code{\link{bridgesde2d}} and \code{\link{bridgesde3d}} functions.  
  \preformatted{
  R> ## 1-dim bridge sde
  R> f <- expression( 2*(1-x) )
  R> g <- expression( 1 )
  R> mod1 <- bridgesde1d(drift=f,diffusion=g,x0=2,y=1,M=1000)
  R> mod1
  Itô Bridges Sde 1D:
        | dX(t) = 2 * (1 - X(t)) * dt + 1 * dW(t)
  Method:
        | Euler scheme with order 0.5
  Summary:
        | Size of process       | N = 1001.
        | Crossing realized     | C = 843 among 1000.
        | Initial value         | x0 = 2.
        | Ending value          | y = 1.
        | Time of process       | t in [0,1].
        | Discretization        | Dt = 0.001.	
  R> summary(mod1) ## Monte-Carlo statistics at T/2=0.5
  
        Monte-Carlo Statistics for X(t) at time t = 0.5
                Crossing realized 843 among 1000                         
   Mean              1.31263
   Variance          0.18352
   Median            1.30504
   Mode              1.46713
   First quartile    1.02722
   Third quartile    1.60984
   Minimum          -0.22080
   Maximum           2.83339
   Skewness          0.01722
   Kurtosis          3.19888
   Coef-variation    0.32636
   3th-order moment  0.00135
   4th-order moment  0.10773
   5th-order moment  0.00645
   6th-order moment  0.11233
  R> plot(mod1)
  R> den <- dsde1d(mod1)
   Density of X(t-t0)|X(t0)=x0 at time t = 1

  Data: x (843 obs.);     Bandwidth 'bw' = 0.2339

             x                f(x)        
  Min.   :0.29822   Min.   :0.01913  
  1st Qu.:0.64911   1st Qu.:0.13600  
  Median :1.00000   Median :0.55258  
  Mean   :1.00000   Mean   :0.70988  
  3rd Qu.:1.35089   3rd Qu.:1.28369  
  Max.   :1.70178   Max.   :1.70511
  R> plot(den)
  
  R> ## 2 and 3-dim Bridge sde
  R> example(bridgesde2d)
  R> example(bridgesde3d) 
  }  
\bold{Estimate the parameters of 1-dim sde:}\cr\cr
Consider a process solution of the general stochastic differential equation:
\deqn{dX_{t} = f(t,X_{t},\underline{\theta}) dt + g(t,X_{t},\underline{\theta}) dW_{t}}{dX(t) = f(t,X(t),theta) dt + g(t,X(t),theta) dW(t)}
The package Sim.DiffProc implements the function \code{\link{fitsde}} of estimate drift and diffusion parameters \eqn{\underline{\theta}=(\theta_{1},\theta_{2},\dots,\theta_{p})}{theta=(theta1,theta2,...,thetap)} 
with different methods of maximum pseudo-likelihood of the 1-dim stochastic differential equation.

An example we use a real data, fit with the CKLS model:
\deqn{dX_{t} = (\theta_{1}+\theta_{2} X_{t}) dt + \theta_{3} X^{\theta_{4}}_{t} dW_{t}}{dX(t) = (theta1+theta2*X(t)) dt + theta3 * X(t)^theta4 dW(t)}
we estimate the vector of parameters \eqn{\underline{\theta}=(\theta_{1},\theta_{2},\theta_{3},\theta_{4})}{theta=(theta1,theta2,theta3,theta4)}, using Euler pseudo-likelihood.
  \preformatted{
  R> ## 1-dim fitsde  
  R> data(Irates)
  R> rates <- Irates[,"r1"]
  R> rates <- window(rates, start=1964.471, end=1989.333)
  R> fx <- expression(theta[1]+theta[2]*x)
  R> gx <- expression(theta[3]*x^theta[4]) 
  R> ## theta = (theta1,theta2,theta3,theta4), p=4
  R> fitmod <- fitsde(rates,drift=fx,diffusion=gx,pmle="euler",start = list(theta1=1,
                      theta2=1,theta3=1,theta4=1),optim.method = "L-BFGS-B")
  R> fitmod
  Call:
  fitsde(data = rates, drift = fx, diffusion = gx, pmle = "euler", 
     start = list(theta1 = 1, theta2 = 1, theta3 = 1, theta4 = 1), 
     optim.method = "L-BFGS-B")
  Coefficients:
      theta1     theta2     theta3     theta4 
   2.0769516 -0.2631871  0.1302158  1.4513173 
  R> summary(fitmod)
  Pseudo maximum likelihood estimation
  Method:  Euler
  Call:
  fitsde(data = rates, drift = fx, diffusion = gx, pmle = "euler", 
     start = list(theta1 = 1, theta2 = 1, theta3 = 1, theta4 = 1), 
     optim.method = "L-BFGS-B")
  Coefficients:
           Estimate Std. Error
  theta1  2.0769516 0.98838467
  theta2 -0.2631871 0.19544290
  theta3  0.1302158 0.02523105
  theta4  1.4513173 0.10323740

  -2 log L: 475.7572 
  R> coef(fitmod)
     theta1     theta2     theta3     theta4 
  2.0769516 -0.2631871  0.1302158  1.4513173 
  R> logLik(fitmod)
  [1] -237.8786
  R> AIC(fitmod)
  [1] 483.7572
  R> BIC(fitmod)
  [1] 487.1514
  R> vcov(fitmod)
                theta1        theta2        theta3        theta4
  theta1  0.9769042534 -1.843596e-01 -2.714334e-04  0.0011374342
  theta2 -0.1843595796  3.819793e-02  5.169849e-05 -0.0002165286
  theta3 -0.0002714334  5.169849e-05  6.366061e-04 -0.0025457493
  theta4  0.0011374342 -2.165286e-04 -2.545749e-03  0.0106579616
  R> confint(fitmod,level=0.95)
              2.5 \%    97.5 \%
  theta1  0.13975321 4.0141499
  theta2 -0.64624812 0.1198740
  theta3  0.08076388 0.1796678
  theta4  1.24897569 1.6536589
  }
\bold{Transition density and Random number generators (RN's) for 1,2 and 3-dim sde:}\cr\cr
Simulation \code{M}-sample for the random variable \eqn{X_{at}}{X(at)} at time \eqn{t=at}{t=at} by a simulated 1, 2 and 3-dim sde, using the functions \code{\link{rsde1d}}, \code{\link{rsde2d}} and \code{\link{rsde3d}}. And \code{\link{dsde1d}}, \code{\link{dsde2d}} and \code{\link{dsde3d}} returns a kernel approximate of transitional densities.   
  \preformatted{
  R> f <- expression(-2*(x<=0)+2*(x>=0))
  R> g <- expression(0.5)
  R> res1 <- snssde1d(drift=f,diffusion=g,M=50,type="str",T=10)
  R> x <- rsde1d(res1,at=10)
  R> x
    [1] -17.64115  21.67111 -19.00162 -20.21546  20.65829  19.59535
    [7] -20.00676 -18.75649 -19.04453 -15.55535 -18.75077  18.89528
   [13] -22.99474 -19.66526 -19.75898  22.02310 -19.68301 -19.08581
   [19] -19.15081 -19.24476 -22.24332  17.74989  19.88449 -18.17091
   [25] -18.65697  19.08473 -17.81218  19.58453  19.27531 -21.88292
   [31]  19.03283 -19.29196  21.99163  20.12123  21.09657 -20.20252
   [37]  20.85097 -19.41987 -18.67530 -19.36289  19.50057  16.30538
   [43]  19.34247 -17.97358  22.81003 -18.40051 -18.47490 -21.86839
   [49] -21.32638 -18.96264
  R> summary(x)
      Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   -22.990 -19.350 -18.440  -4.436  19.330  22.810
  R> den <- dsde1d(res1,at=10)
  R> ## 2 and 3-dim rsde
  R> example(dsde2d)
  R> example(dsde3d)
  }
\bold{First-passage-time (f.p.t) in 1,2 and 3-dim sde}\cr\cr
The functions \code{\link{fptsde1d}} (\code{\link{fptsde2d}} and \code{\link{fptsde3d}} for 2 and 3-dim) returns a random variable  \eqn{\tau_{(X(t),S(t))}}{tau(X(t),S(t))} "first passage time", is defined as: 
\deqn{\tau_{(X(t),S(t))} = \{ t \geq 0 ; X_{t} \geq S(t) \},\quad if \quad X(t_{0}) < S(t_{0})}{tau(X(t),S(t))={t>=0; X(t) >= S(t)}, if X(t0) < S(t0)}
\deqn{\tau_{(X(t),S(t))} = \{ t \geq 0 ; X_{t} \leq S(t) \},\quad if \quad X(t_{0}) > S(t_{0})}{tau(X(t),S(t))={t>=0; X(t) <= S(t)}, if X(t0) > S(t0)}
And \code{\link{dfptsde1d}}, \code{\link{dfptsde2d}} and \code{\link{dfptsde3d}} returns a kernel density approximation for first passage time.
with \eqn{S(t)}{S(t)} is through a continuous boundary (barrier). 
  \preformatted{
  R> f <- expression( 0.5*x*t )
  R> g <- expression( sqrt(1+x^2) )
  R> St <- expression(-0.5*sqrt(t)+exp(t^2))
  R> mod <- snssde1d(drift=f,diffusion=g,x0=2,M=1000)
  R> fptmod <- fptsde1d(mod,boundary=St)
  R> fptmod
   Itô Sde 1D:
        | dX(t) = 0.5 * X(t) * t * dt + sqrt(1 + X(t)^2) * dW(t)
        | t in [0,1].
   Boundary:
        | S(t) = -0.5 * sqrt(t) + exp(t^2)
   F.P.T:
        | T(S(t),X(t)) = inf{t >=  0 : X(t) <=  -0.5 * sqrt(t) + exp(t^2) }
        | Crossing realized 738 among 1000.
  R> summary(fptmod)

   Monte-Carlo Statistics of F.P.T:
   |T(S(t),X(t)) = inf{t >=  0 : X(t) <=  -0.5 * sqrt(t) + exp(t^2) }
                        
   Mean             0.47742
   Variance         0.07348
   Median           0.44831
   Mode             0.18582
   First quartile   0.23746
   Third quartile   0.71321
   Minimum          0.03002
   Maximum          0.98877
   Skewness         0.22793
   Kurtosis         1.79959
   Coef-variation   0.56778
   3th-order moment 0.00454
   4th-order moment 0.00972
   5th-order moment 0.00134
   6th-order moment 0.00158
  R> den <- dfptsde1d(mod,boundary=St)
  R> den
    Kernel density for the F.P.T of X(t)
    T(S,X) = inf{t >= 0 : X(t) <= -0.5 * sqrt(t) + exp(t^2)}

     Data: fpt (738 obs.);    Bandwidth 'bw' = 0.0828

          x                f(x)       
   Min.   :-0.2095   Min.   :0.0019  
   1st Qu.: 0.1458   1st Qu.:0.2163  
   Median : 0.5010   Median :0.5307  
   Mean   : 0.5010   Mean   :0.7029  
   3rd Qu.: 0.8563   3rd Qu.:1.1943  
   Max.   : 1.2116   Max.   :1.8548
  R> ## fpt in 2 and 3-dim sde
  R> example(dfptsde2d)
  R> example(dfptsde3d)  
  }
  
For other examples see \code{demo(Sim.DiffProc)}, and for an overview of this package, see \code{\link[utils]{browseVignettes}('Sim.DiffProc')} for more informations. 
}
\section{Requirements}{
\R version >= 3.0.0
}
\author{
A.C. Guidoum \email{acguidoum@usthb.dz} and K. Boukhetala \email{kboukhetala@usthb.dz} (Dept. Probability and Statistics, USTHB, Algeria).\cr 
Please send comments, error reports, etc. to the author via the addresses email.
}
\section{Licence}{
This package and its documentation are usable under the terms of the "GNU
General Public License", a copy of which is distributed with the package.
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

   SaItô, Y, and Mitsui, T. (1993).
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
sde, yumia, QPot, DiffusionRgqd, fptdApprox.
}

