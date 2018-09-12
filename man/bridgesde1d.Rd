

\name{bridgesde1d}
\alias{bridgesde1d}
\alias{bridgesde1d.default}
\alias{print.bridgesde1d}
\alias{time.bridgesde1d}
\alias{mean.bridgesde1d}
\alias{summary.bridgesde1d}
\alias{Median.bridgesde1d}
\alias{Mode.bridgesde1d}
\alias{quantile.bridgesde1d}
\alias{kurtosis.bridgesde1d}
\alias{skewness.bridgesde1d}
\alias{moment.bridgesde1d}
\alias{cv.bridgesde1d}
\alias{max.bridgesde1d}
\alias{min.bridgesde1d}
\alias{bconfint.bridgesde1d}
\alias{plot.bridgesde1d}
\alias{points.bridgesde1d}
\alias{lines.bridgesde1d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulation of 1-D Bridge SDE}
\description{The (S3) generic function \code{bridgesde1d} for simulation of 1-dim bridge stochastic differential equations,Itô or Stratonovich type, with different methods.}
\usage{
bridgesde1d(N, \dots)
\method{bridgesde1d}{default}(N = 1000, M=1, x0 = 0, y = 0, t0 = 0, T = 1, Dt=NULL, 
   drift, diffusion, alpha = 0.5, mu = 0.5, type = c("ito", "str"), 
   method = c("euler", "milstein", "predcorr", "smilstein", "taylor", 
   "heun", "rk1", "rk2", "rk3"), \dots)
   
\method{summary}{bridgesde1d}(object, at ,digits=NULL, \dots)							  
\method{time}{bridgesde1d}(x, \dots)
\method{mean}{bridgesde1d}(x, at, \dots)
\method{Median}{bridgesde1d}(x, at, \dots)
\method{Mode}{bridgesde1d}(x, at, \dots)
\method{quantile}{bridgesde1d}(x, at, \dots)
\method{kurtosis}{bridgesde1d}(x, at, \dots)
\method{skewness}{bridgesde1d}(x, at, \dots)
\method{min}{bridgesde1d}(x, at, \dots)
\method{max}{bridgesde1d}(x, at, \dots)
\method{moment}{bridgesde1d}(x, at, \dots)
\method{cv}{bridgesde1d}(x, at, \dots)
\method{bconfint}{bridgesde1d}(x, at,  \dots)

\method{plot}{bridgesde1d}(x, \dots)
\method{lines}{bridgesde1d}(x, \dots)
\method{points}{bridgesde1d}(x, \dots)							  
}
%- maybe also 'usage' for other xs documented here.
\arguments{
  \item{N}{number of simulation steps.}
  \item{M}{number of trajectories.}
  \item{x0}{initial value of the process at time \code{t0}.}
  \item{y}{terminal value of the process at time \code{T}.}
  \item{t0}{initial time.}
  \item{T}{final time.} 
  \item{Dt}{time step of the simulation (discretization). If it is \code{\link[base]{NULL}} a default \eqn{\Delta t = \frac{T-t_{0}}{N}}{Dt = (T-t0)/N}.}
  \item{drift}{drift coefficient: an \code{\link[base]{expression}} of two variables \code{t} and \code{x}.}
  \item{diffusion}{diffusion coefficient: an \code{\link[base]{expression}} of two variables \code{t} and \code{x}.}
  \item{alpha, mu}{weight of the predictor-corrector scheme; the default \code{alpha = 0.5} and \code{mu = 0.5}.}
  \item{type}{if \code{type="ito"} simulation diffusion bridge of Itô type, else \code{type="str"} simulation diffusion bridge of Stratonovich type; the default \code{type="ito"}.}
  \item{method}{numerical methods of simulation, the default \code{method = "euler"}; see \code{\link{snssde1d}}.}
  \item{x, object}{an object inheriting from class \code{"bridgesde1d"}.}
  \item{at}{time between \code{t0} and \code{T}. Monte-Carlo statistics of the solution \eqn{X_{t}}{X(t)} at time \code{at}. The default \code{at = T/2}.} 
  \item{digits}{integer, used for number formatting.}  
  \item{\dots}{potentially further arguments for (non-default) methods.}
}
\details{
The function \code{bridgesde1d} returns a trajectory of the diffusion bridge starting at \code{x} at time \code{t0} and ending at \code{y} at time \code{T}.

The methods of approximation are classified according to their different properties. Mainly two criteria of optimality are used in the literature: the strong 
and the weak (orders of) convergence. The \code{method} of simulation can be one among: \code{Euler-Maruyama Order 0.5}, \code{Milstein Order 1}, \code{Milstein Second-Order}, 
\code{Predictor-Corrector method}, \code{Itô-Taylor Order 1.5}, \code{Heun Order 2} and \code{Runge-Kutta Order 1, 2 and 3}.

An overview of this package, see \code{\link[utils]{browseVignettes}('Sim.DiffProc')} for more informations. 
}
\value{
\code{bridgesde1d} returns an object inheriting from \code{\link[base]{class}} \code{"bridgesde1d"}.
\item{X }{an invisible \code{\link[stats]{ts}} object.}
\item{drift }{drift coefficient.}
\item{diffusion }{diffusion coefficient.}
\item{C}{indices of crossing realized of X(t).}
\item{type }{type of sde.}
\item{method }{the numerical method used.}
}
\references{
   Bladt, M. and Sorensen, M. (2007).
   Simple simulation of diffusion bridges with application to likelihood inference for diffusions.
   \emph{Working Paper, University of Copenhagen}. 

   Iacus, S.M. (2008). 
   \emph{Simulation and inference for stochastic differential equations: with R examples}. 
   Springer-Verlag, New York
}
\author{
A.C. Guidoum, K. Boukhetala.
}
\seealso{
\code{\link{bridgesde2d}} and \code{\link{bridgesde3d}} for 2 and 3-dim. 

\code{\link[sde]{DBridge}} in package "sde". 
}
\examples{
## Example 1: Ito bridge sde
## Ito Bridge sde
## dX(t) = 2*(1-X(t)) *dt + dW(t)
## x0 = 2 at time t0=0 , and y = 1 at time T=1
set.seed(1234)

f <- expression( 2*(1-x) )
g <- expression( 1 )
mod1 <- bridgesde1d(drift=f,diffusion=g,x0=2,y=1,M=1000)
mod1
summary(mod1) ## Monte-Carlo statistics at T/2=0.5
summary(mod1,at=0.75) ## Monte-Carlo statistics at 0.75
\dontrun{
plot(mod1)
lines(time(mod1),apply(mod1$X,1,mean),col=2,lwd=2)
lines(time(mod1),apply(mod1$X,1,bconfint,level=0.95)[1,],col=4,lwd=2)
lines(time(mod1),apply(mod1$X,1,bconfint,level=0.95)[2,],col=4,lwd=2)
legend("topleft",c("mean path",paste("bound of", 95," percent confidence")),
       inset = .01,col=c(2,4),lwd=2,cex=0.8)
}

## Example 2: Stratonovich sde
## dX(t) = ((2-X(t))/(2-t)) dt + X(t) o dW(t)
## x0 = 2 at time t0=0 , and y = 2 at time T=1
set.seed(1234)

f <- expression((2-x)/(2-t))
g <- expression(x)
mod2 <- bridgesde1d(type="str",drift=f,diffusion=g,M=1000,x0=2,y=2)
mod2
summary(mod2,at = 0.25) ## Monte-Carlo statistics at 0.25
summary(mod2,at = 0.5)  ## Monte-Carlo statistics at 0.5
summary(mod2,at = 0.75 )## Monte-Carlo statistics at 0.75
\dontrun{
plot(mod2)
lines(time(mod2),apply(mod2$X,1,mean),col=2,lwd=2)
lines(time(mod2),apply(mod2$X,1,bconfint,level=0.95)[1,],col=4,lwd=2)
lines(time(mod2),apply(mod2$X,1,bconfint,level=0.95)[2,],col=4,lwd=2)
legend("topright",c("mean path",paste("bound of", 95," percent confidence")),
       inset = .01,col=c(2,4),lwd=2,cex=0.8)
}
}

\keyword{sde}
\keyword{ts}
\keyword{mts}
