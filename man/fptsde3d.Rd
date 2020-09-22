

\name{fptsde3d}
\alias{fptsde3d}
\alias{fptsde3d.default}
\alias{summary.fptsde3d}
\alias{print.fptsde3d}
\alias{mean.fptsde3d}
\alias{Median.fptsde3d}
\alias{Mode.fptsde3d}
\alias{quantile.fptsde3d}
\alias{kurtosis.fptsde3d}
\alias{skewness.fptsde3d}
\alias{moment.fptsde3d}
\alias{cv.fptsde3d}
\alias{max.fptsde3d}
\alias{min.fptsde3d}
\alias{dfptsde3d}
\alias{dfptsde3d.default}
\alias{print.dfptsde3d}
\alias{plot.dfptsde3d}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Approximate densities and random generation for first passage time in 3-D SDE's}
\description{Kernel density and random generation for first-passage-time (f.p.t) in 3-dim stochastic differential equations.}
\usage{
fptsde3d(object, \dots)
dfptsde3d(object, \dots)

\method{fptsde3d}{default}(object, boundary, \dots)
\method{summary}{fptsde3d}(object, digits=NULL, \dots)
\method{mean}{fptsde3d}(x, \dots)
\method{Median}{fptsde3d}(x, \dots)
\method{Mode}{fptsde3d}(x, \dots)
\method{quantile}{fptsde3d}(x, \dots)
\method{kurtosis}{fptsde3d}(x, \dots)
\method{skewness}{fptsde3d}(x, \dots)
\method{min}{fptsde3d}(x, \dots)
\method{max}{fptsde3d}(x, \dots)
\method{moment}{fptsde3d}(x, \dots)
\method{cv}{fptsde3d}(x, \dots)

\method{dfptsde3d}{default}(object, pdf=c("Joint","Marginal"), \dots)
\method{plot}{dfptsde3d}(x,display="rgl",hist=FALSE, \dots)						   
}
%- maybe also 'usage' for other xs documented here.
\arguments{
  \item{object}{an object inheriting from class \code{\link{snssde3d}} for \code{fptsde3d}, and \code{\link{fptsde3d}} for \code{dfptsde3d}.}
  \item{boundary}{an \code{\link[base]{expression}} of a constant or time-dependent boundary.}
  \item{pdf}{probability density function \code{Joint} or \code{Marginal}.} 
  \item{x}{an object inheriting from class \code{dfptsde3d}.}
  \item{digits}{integer, used for number formatting.}
  \item{display}{display plots.} 
  \item{hist}{if \code{hist=TRUE} plot histogram. Based on \code{\link[MASS]{truehist}} function.} 
  \item{\dots}{potentially arguments to be passed to methods, such as \code{\link[stats]{density}} for marginal density  and \code{\link[sm]{sm.density}} for joint density.}
}
\details{
The function \code{fptsde3d} returns a random variable \eqn{(\tau_{(X(t),S(t))},\tau_{(Y(t),S(t))},\tau_{(Z(t),S(t))})}{(tau(X(t),S(t)),tau(Y(t),S(t)),tau(Z(t),S(t)))} "first passage time", is defined as : 
\deqn{\tau_{(X(t),S(t))} = \{ t \geq 0 ; X_{t} \geq S(t) \},\quad if \quad X(t_{0}) < S(t_{0})}{tau(X(t),S(t))={t>=0; X(t) >= S(t)}, if X(t0) < S(t0)}
\deqn{\tau_{(Y(t),S(t))} = \{ t \geq 0 ; Y_{t} \geq S(t) \},\quad if \quad Y(t_{0}) < S(t_{0})}{tau(Y(t),S(t))={t>=0; Y(t) >= S(t)}, if Y(t0) < S(t0)}
\deqn{\tau_{(Z(t),S(t))} = \{ t \geq 0 ; Z_{t} \geq S(t) \},\quad if \quad Z(t_{0}) < S(t_{0})}{tau(Z(t),S(t))={t>=0; Z(t) >= S(t)}, if Z(t0) < S(t0)}
and:
\deqn{\tau_{(X(t),S(t))} = \{ t \geq 0 ; X_{t} \leq S(t) \},\quad if \quad X(t_{0}) > S(t_{0})}{tau(X(t),S(t))={t>=0; X(t) <= S(t)}, if X(t0) > S(t0)}
\deqn{\tau_{(Y(t),S(t))} = \{ t \geq 0 ; Y_{t} \leq S(t) \},\quad if \quad Y(t_{0}) > S(t_{0})}{tau(Y(t),S(t))={t>=0; Y(t) <= S(t)}, if Y(t0) > S(t0)}
\deqn{\tau_{(Z(t),S(t))} = \{ t \geq 0 ; Z_{t} \leq S(t) \},\quad if \quad Z(t_{0}) > S(t_{0})}{tau(Z(t),S(t))={t>=0; Z(t) <= S(t)}, if Z(t0) > S(t0)}

\if{html}{\figure{fpt3d.png}{fig11}}

And \code{dfptsde3d} returns a marginal kernel density approximation for \eqn{(\tau_{(X(t),S(t))},\tau_{(Y(t),S(t))},\tau_{(Z(t),S(t))})}{(tau(X(t),S(t)),tau(Y(t),S(t)),tau(Z(t),S(t)))} "first passage time".
with \eqn{S(t)}{S(t)} is through a continuous boundary (barrier).

\if{html}{\figure{dfpt3d.png}{fig12}}

An overview of this package, see \code{\link[utils]{browseVignettes}('Sim.DiffProc')} for more informations.
}
\value{
\code{dfptsde3d} {gives the marginal kernel density approximation for fpt.}
\code{fptsde3d} {generates random of fpt.}
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
\code{\link{fptsde1d}} for simulation fpt in sde 1-dim. \code{\link{fptsde2d}} for simulation fpt in sde 2-dim.

\code{FPTL} for computes values of the first passage time location (FPTL) function, and \code{Approx.fpt.density} 
for approximate first-passage-time (f.p.t.) density in package "fptdApprox". 

\code{GQD.TIpassage} for compute the First Passage Time Density of a GQD With Time Inhomogeneous Coefficients in package "DiffusionRgqd".

}
\examples{

## dX(t) = 4*(-1-X(t))*Y(t) dt + 0.2 * dW1(t) 
## dY(t) = 4*(1-Y(t)) *X(t) dt + 0.2 * dW2(t) 
## dZ(t) = 4*(1-Z(t)) *Y(t) dt + 0.2 * dW3(t) 
## x0 = 0, y0 = -2, z0 = 0, and barrier -3+5*t.       
## W1(t), W2(t) and W3(t) three independent Brownian motion      
set.seed(1234)

# SDE's 3d

fx <- expression(4*(-1-x)*y, 4*(1-y)*x, 4*(1-z)*y)
gx <- rep(expression(0.2),3)
mod3d <- snssde3d(drift=fx,diffusion=gx,M=500)

# boundary 
St <- expression(-3+5*t)

# random

out <- fptsde3d(mod3d,boundary=St)
out
summary(out)

# Marginal density

denM <- dfptsde3d(out,pdf="M")
denM
plot(denM)

# Multiple isosurfaces
\dontrun{
denJ <- dfptsde3d(out,pdf="J")
denJ
plot(denJ,display="rgl")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{fpt}
\keyword{sde}
\keyword{ts}
\keyword{mts}
