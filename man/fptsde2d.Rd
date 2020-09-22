

\name{fptsde2d}
\alias{fptsde2d}
\alias{fptsde2d.default}
\alias{summary.fptsde2d}
\alias{print.fptsde2d}
\alias{mean.fptsde2d}
\alias{Median.fptsde2d}
\alias{Mode.fptsde2d}
\alias{quantile.fptsde2d}
\alias{kurtosis.fptsde2d}
\alias{skewness.fptsde2d}
\alias{moment.fptsde2d}
\alias{cv.fptsde2d}
\alias{max.fptsde2d}
\alias{min.fptsde2d}
\alias{dfptsde2d}
\alias{dfptsde2d.default}
\alias{print.dfptsde2d}
\alias{plot.dfptsde2d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Approximate densities and random generation for first passage time in 2-D SDE's}
\description{Kernel density and random generation for first-passage-time (f.p.t) in 2-dim stochastic differential equations.}
\usage{
fptsde2d(object, \dots)
dfptsde2d(object, \dots)

\method{fptsde2d}{default}(object, boundary, \dots)
\method{summary}{fptsde2d}(object, digits=NULL, \dots)
\method{mean}{fptsde2d}(x, \dots)
\method{Median}{fptsde2d}(x, \dots)
\method{Mode}{fptsde2d}(x, \dots)
\method{quantile}{fptsde2d}(x, \dots)
\method{kurtosis}{fptsde2d}(x, \dots)
\method{skewness}{fptsde2d}(x, \dots)
\method{min}{fptsde2d}(x, \dots)
\method{max}{fptsde2d}(x, \dots)
\method{moment}{fptsde2d}(x, \dots)
\method{cv}{fptsde2d}(x, \dots)

\method{dfptsde2d}{default}(object, pdf=c("Joint","Marginal"), \dots)
\method{plot}{dfptsde2d}(x,display=c("persp","rgl","image","contour"),
                         hist=FALSE, \dots)						  
}
%- maybe also 'usage' for other xs documented here.
\arguments{
  \item{object}{an object inheriting from class \code{\link{snssde2d}} for \code{fptsde2d}, and \code{\link{fptsde2d}} for \code{dfptsde2d}.}
  \item{boundary}{an \code{\link[base]{expression}} of a constant or time-dependent boundary.} 
  \item{pdf}{probability density function \code{Joint} or \code{Marginal}.}    
  \item{x}{an object inheriting from class \code{fptsde2d}.}
  \item{digits}{integer, used for number formatting.}
  \item{display}{display plots.} 
  \item{hist}{if \code{hist=TRUE} plot histogram. Based on \code{\link[MASS]{truehist}} function.} 
  \item{\dots}{potentially further arguments for (non-default) methods. arguments to be passed to methods, such as \code{\link[stats]{density}} for marginal density and \code{\link[MASS]{kde2d}} fro joint density.}
}
\details{
The function \code{fptsde1d} returns a random variable \eqn{(\tau_{(X(t),S(t))},\tau_{(Y(t),S(t))})}{(tau(X(t),S(t)),tau(Y(t),S(t)))} "first passage time", is defined as : 
\deqn{\tau_{(X(t),S(t))} = \{ t \geq 0 ; X_{t} \geq S(t) \},\quad if \quad X(t_{0}) < S(t_{0})}{tau(X(t),S(t))={t>=0; X(t) >= S(t)}, if X(t0) < S(t0)}
\deqn{\tau_{(Y(t),S(t))} = \{ t \geq 0 ; Y_{t} \geq S(t) \},\quad if \quad Y(t_{0}) < S(t_{0})}{tau(Y(t),S(t))={t>=0; Y(t) >= S(t)}, if Y(t0) < S(t0)}
and:
\deqn{\tau_{(X(t),S(t))} = \{ t \geq 0 ; X_{t} \leq S(t) \},\quad if \quad X(t_{0}) > S(t_{0})}{tau(X(t),S(t))={t>=0; X(t) <= S(t)}, if X(t0) > S(t0)}
\deqn{\tau_{(Y(t),S(t))} = \{ t \geq 0 ; Y_{t} \leq S(t) \},\quad if \quad Y(t_{0}) > S(t_{0})}{tau(Y(t),S(t))={t>=0; Y(t) <= S(t)}, if Y(t0) > S(t0)}

\if{html}{\figure{fpt2d.png}{fig09}}

And \code{dfptsde2d} returns a kernel density approximation for \eqn{(\tau_{(X(t),S(t))},\tau_{(Y(t),S(t))})}{(tau(X(t),S(t)),tau(Y(t),S(t)))} "first passage time".
with \eqn{S(t)}{S(t)} is through a continuous boundary (barrier).

\if{html}{\figure{dfpt2d.png}{fig10}}

An overview of this package, see \code{\link[utils]{browseVignettes}('Sim.DiffProc')} for more informations.
}
\value{
\code{dfptsde2d} {gives the kernel density approximation for fpt.}
\code{fptsde2d} {generates random of fpt.}
}
\references{
  Argyrakisa, P. and G.H. Weiss (2006).
  A first-passage time problem for many random walkers.
  \emph{Physica A}. \bold{363}, 343--347.
  
  Aytug H., G. J. Koehler (2000). 
  New stopping criterion for genetic algorithms. 
  \emph{European Journal of Operational Research}, \bold{126}, 662--674.
  
  Boukhetala, K. (1996)
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
  
   Roman, R.P., Serrano, J. J., Torres, F. (2008).
   First-passage-time location function: Application to determine first-passage-time densities in diffusion processes.
   \emph{Computational Statistics and Data Analysis}. \bold{52}, 4132--4146.
  
   Roman, R.P., Serrano, J. J., Torres, F. (2012). 
   An R package for an efficient approximation of first-passage-time densities for diffusion processes based on the FPTL function. 
   \emph{Applied Mathematics and Computation}, \bold{218}, 8408--8428.
    
  Gardiner, C. W. (1997). 
  \emph{Handbook of Stochastic Methods}.
  Springer-Verlag, New York.
}
\author{
A.C. Guidoum, K. Boukhetala.
}

\seealso{
\code{\link[Sim.DiffProc]{fptsde1d}} for simulation fpt in sde 1-dim. \code{\link{fptsde3d}} for simulation fpt in sde 3-dim. 

\code{FPTL} for computes values of the first passage time location (FPTL) function, and \code{Approx.fpt.density} 
for approximate first-passage-time (f.p.t.) density in package "fptdApprox". 

\code{GQD.TIpassage} for compute the First Passage Time Density of a GQD With Time Inhomogeneous Coefficients in package "DiffusionRgqd".

}
\examples{

## dX(t) = 5*(-1-Y(t))*X(t) * dt + 0.5 * dW1(t)          
## dY(t) = 5*(-1-X(t))*Y(t) * dt + 0.5 * dW2(t)
## x0 = 2, y0 = -2, and barrier -3+5*t.
## W1(t) and W2(t) two independent Brownian motion
set.seed(1234)

# SDE's 2d
fx <- expression(5*(-1-y)*x , 5*(-1-x)*y)
gx <- expression(0.5 , 0.5)
mod2d <- snssde2d(drift=fx,diffusion=gx,x0=c(2,-2),M=100)

# boundary

St <- expression(-1+5*t)

# random fpt

out <- fptsde2d(mod2d,boundary=St)
out
summary(out)

# Marginal density 

denM <- dfptsde2d(out,pdf="M")
denM
plot(denM)

# Joint density

denJ <- dfptsde2d(out,pdf="J",n=200,lims=c(0.28,0.4,0.04,0.13))
denJ
plot(denJ)
plot(denJ,display="image")
plot(denJ,display="image",drawpoints=TRUE,cex=0.5,pch=19,col.pt='green')
plot(denJ,display="contour")
plot(denJ,display="contour",color.palette=colorRampPalette(c('white','green','blue','red')))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{fpt}
\keyword{sde}
\keyword{ts}
\keyword{mts}
