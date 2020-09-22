

\name{fptsde1d}
\alias{fptsde1d}
\alias{fptsde1d.default}
\alias{summary.fptsde1d}
\alias{print.fptsde1d}
\alias{mean.fptsde1d}
\alias{Median.fptsde1d}
\alias{Mode.fptsde1d}
\alias{quantile.fptsde1d}
\alias{kurtosis.fptsde1d}
\alias{skewness.fptsde1d}
\alias{moment.fptsde1d}
\alias{cv.fptsde1d}
\alias{max.fptsde1d}
\alias{min.fptsde1d}
\alias{dfptsde1d}
\alias{dfptsde1d.default}
\alias{print.dfptsde1d}
\alias{plot.dfptsde1d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Approximate densities and random generation for first passage time in 1-D SDE}
\description{Kernel density and random generation for first-passage-time (f.p.t) in 1-dim stochastic differential equations.}
\usage{
fptsde1d(object, \dots)
dfptsde1d(object, \dots)

\method{fptsde1d}{default}(object, boundary, \dots)
\method{summary}{fptsde1d}(object, digits=NULL, \dots)
\method{mean}{fptsde1d}(x, \dots)
\method{Median}{fptsde1d}(x, \dots)
\method{Mode}{fptsde1d}(x, \dots)
\method{quantile}{fptsde1d}(x, \dots)
\method{kurtosis}{fptsde1d}(x, \dots)
\method{skewness}{fptsde1d}(x, \dots)
\method{min}{fptsde1d}(x, \dots)
\method{max}{fptsde1d}(x, \dots)
\method{moment}{fptsde1d}(x, \dots)
\method{cv}{fptsde1d}(x, \dots)

\method{dfptsde1d}{default}(object, \dots)
\method{plot}{dfptsde1d}(x, hist=FALSE, \dots)						  
}
%- maybe also 'usage' for other xs documented here.
\arguments{
  \item{object}{an object inheriting from class \code{\link{snssde1d}} for \code{fptsde1d}, and \code{\link{fptsde1d}} for \code{dfptsde1d}.}
  \item{boundary}{an \code{\link[base]{expression}} of a constant or time-dependent boundary.}  
  \item{x}{an object inheriting from class \code{dfptsde1d}.}
  \item{hist}{if \code{hist=TRUE} plot histogram. Based on \code{\link[MASS]{truehist}} function.} 
  \item{digits}{integer, used for number formatting.}  
  \item{\dots}{potentially further arguments for (non-default) methods, such as \code{\link[stats]{density}} for \code{dfptsde1d}.}
}
\details{
The function \code{fptsde1d} returns a random variable \eqn{\tau_{(X(t),S(t))}}{tau(X(t),S(t))} "first passage time", is defined as : 
\deqn{\tau_{(X(t),S(t))} = \{ t \geq 0 ; X_{t} \geq S(t) \},\quad if \quad X(t_{0}) < S(t_{0})}{tau(X(t),S(t))={t>=0; X(t) >= S(t)}, if X(t0) < S(t0)}
\deqn{\tau_{(X(t),S(t))} = \{ t \geq 0 ; X_{t} \leq S(t) \},\quad if \quad X(t_{0}) > S(t_{0})}{tau(X(t),S(t))={t>=0; X(t) <= S(t)}, if X(t0) > S(t0)}

\if{html}{\figure{fpt1d.png}{fig07}}

And \code{dfptsde1d} returns a kernel density approximation for \eqn{\tau_{(X(t),S(t))}}{tau(X(t),S(t))} "first passage time".
with \eqn{S(t)}{S(t)} is through a continuous boundary (barrier).

\if{html}{\figure{dfpt1d.png}{fig08}}

An overview of this package, see \code{\link[utils]{browseVignettes}('Sim.DiffProc')} for more informations. 
}
\value{
\code{dfptsde1d} {gives the density estimate of fpt.}
\code{fptsde1d} {generates random of fpt.}
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
\code{\link[Sim.DiffProc]{fptsde2d}} and \code{\link[Sim.DiffProc]{fptsde3d}} simulation fpt for 2 and 3-dim SDE.

\code{FPTL} for computes values of the first passage time location (FPTL) function, and \code{Approx.fpt.density} 
for approximate first-passage-time (f.p.t.) density in package "fptdApprox". 

\code{GQD.TIpassage} for compute the First Passage Time Density of a GQD With Time Inhomogeneous Coefficients in package "DiffusionRgqd".

}
\examples{

## Example 1: Ito SDE
## dX(t) = -4*X(t) *dt + 0.5*dW(t)
## S(t) = 0 (constant boundary)
set.seed(1234)

# SDE 1d
f <- expression( -4*x )
g <- expression( 0.5 )
mod <- snssde1d(drift=f,diffusion=g,x0=2,M=1000)

# boundary
St <- expression(0) 

# random
out <- fptsde1d(mod, boundary=St)
out
summary(out)
# density approximate
den <- dfptsde1d(out)
den
plot(den)

## Example 2: Stratonovich SDE
## dX(t) = 0.5*X(t)*t *dt + sqrt(1+X(t)^2) o dW(t)
## S(t) = -0.5*sqrt(t) + exp(t^2) (time-dependent boundary)
set.seed(1234)

# SDE 1d
f <- expression( 0.5*x*t )
g <- expression( sqrt(1+x^2) )
mod2 <- snssde1d(drift=f,diffusion=g,x0=2,M=1000,type="srt")

# boundary
St <- expression(-0.5*sqrt(t)+exp(t^2))

# random
out2 <- fptsde1d(mod2,boundary=St)
out2
summary(out2)
# density approximate
plot(dfptsde1d(out2,bw='ucv'))


## Example 3: fptsde1d vs fptdApproximate
\dontrun{
f <- expression( -0.5*x+0.5*5 )
g <- expression( 1 )
St <- expression(5+0.25*sin(2*pi*t)) 
mod <- snssde1d(drift=f,diffusion=g,boundary=St,x0=3,T=10,N=10^4,M =10000)
mod

# random
out3 <- fptsde1d(mod,boundary=St)
out3
summary(out3)
# density approximate: 
library("fptdApprox")
# Under `fptdApprox':
# Define the diffusion process and give its transitional density:
OU <- diffproc(c("alpha*x + beta","sigma^2",
"dnorm((x-(y*exp(alpha*(t-s)) - beta*(1 - exp(alpha*(t-s)))/alpha))/
(sigma*sqrt((exp(2*alpha*(t-s)) - 1)/(2*alpha))),0,1)/
(sigma*sqrt((exp(2*alpha*(t-s)) - 1)/(2*alpha)))",
"pnorm(x, y*exp(alpha*(t-s)) - beta*(1 - exp(alpha*(t-s)))/alpha,
sigma*sqrt((exp(2*alpha*(t-s)) - 1)/(2*alpha)))"))
# Approximate the first passgage time density for OU, starting in X_0 = 3
# passing through 5+0.25*sin(2*pi*t) on the time interval [0,10]:
res <- Approx.fpt.density(OU, 0, 10, 3,"5+0.25*sin(2*pi*t)", list(alpha=-0.5,beta=0.5*5,sigma=1))

##
plot(dfptsde1d(out3,bw='ucv'),main = 'fptsde1d vs fptdApproximate')
lines(res$y~res$x, type = 'l',lwd=2)
legend('topright', lty = c('solid', 'dashed'), col = c(1, 2),
       legend = c('fptdApproximate', 'fptsde1d'), lwd = 2, bty = 'n')

}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{fpt}
\keyword{sde}
\keyword{ts}
