

\name{bridgesde2d}
\alias{bridgesde2d}
\alias{bridgesde2d.default}
\alias{print.bridgesde2d}
\alias{time.bridgesde2d}
\alias{mean.bridgesde2d}
\alias{Median.bridgesde2d}
\alias{Mode.bridgesde2d}
\alias{quantile.bridgesde2d}
\alias{kurtosis.bridgesde2d}
\alias{skewness.bridgesde2d}
\alias{moment.bridgesde2d}
\alias{summary.bridgesde2d}
\alias{cv.bridgesde2d}
\alias{max.bridgesde2d}
\alias{min.bridgesde2d}
\alias{bconfint.bridgesde2d}
\alias{plot.bridgesde2d}
\alias{points.bridgesde2d}
\alias{lines.bridgesde2d}
\alias{plot2d.bridgesde2d}
\alias{points2d.bridgesde2d}
\alias{lines2d.bridgesde2d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulation of 2-D Bridge SDE's}
\description{The (S3) generic function \code{bridgesde2d} for simulation of 2-dim bridge stochastic differential equations,Itô or Stratonovich type, with different methods.}
\usage{
bridgesde2d(N, \dots)
\method{bridgesde2d}{default}(N = 1000, M = 1, x0 = c(0, 0), 
   y = c(0, 0),t0 = 0, T = 1, Dt,drift, diffusion, corr = NULL,
   alpha = 0.5, mu = 0.5,type = c("ito", "str"),method = 
   c("euler", "milstein","predcorr", "smilstein", "taylor",
   "heun", "rk1", "rk2", "rk3"), \dots)

\method{summary}{bridgesde2d}(object, at,
     digits=NULL, \dots)								  
\method{time}{bridgesde2d}(x, \dots)
\method{mean}{bridgesde2d}(x, at, \dots)
\method{Median}{bridgesde2d}(x, at, \dots)
\method{Mode}{bridgesde2d}(x, at, \dots)
\method{quantile}{bridgesde2d}(x, at, \dots)
\method{kurtosis}{bridgesde2d}(x, at, \dots)
\method{skewness}{bridgesde2d}(x, at, \dots)
\method{min}{bridgesde2d}(x, at, \dots)
\method{max}{bridgesde2d}(x, at, \dots)
\method{moment}{bridgesde2d}(x, at, \dots)
\method{cv}{bridgesde2d}(x, at, \dots)
\method{bconfint}{bridgesde2d}(x, at, \dots)

\method{plot}{bridgesde2d}(x, \dots)
\method{lines}{bridgesde2d}(x, \dots)
\method{points}{bridgesde2d}(x, \dots)	
\method{plot2d}{bridgesde2d}(x, \dots)
\method{lines2d}{bridgesde2d}(x, \dots)
\method{points2d}{bridgesde2d}(x, \dots)								  
}
%- maybe also 'usage' for other xs documented here.
\arguments{
  \item{N}{number of simulation steps.}
  \item{M}{number of trajectories.}
  \item{x0}{initial value (numeric vector of length 2) of the process \eqn{X_t}{X(t)} and \eqn{Y_t}{Y(t)} at time \eqn{t_0}{t0}.}
  \item{y}{terminal value (numeric vector of length 2) of the process \eqn{X_t}{X(t)} and \eqn{Y_t}{Y(t)} at time \eqn{T}{T}.}
  \item{t0}{initial time.}
  \item{T}{final time.} 
  \item{Dt}{time step of the simulation (discretization). If it is \code{\link[base]{missing}} a default \eqn{\Delta t = \frac{T-t_{0}}{N}}{Dt = (T-t0)/N}.}
  \item{drift}{drift coefficient: an \code{\link[base]{expression}} of three variables \code{t}, \code{x} and \code{y} for process \eqn{X_t}{X(t)} and \eqn{Y_t}{Y(t)}.}
  \item{diffusion}{diffusion coefficient: an \code{\link[base]{expression}} of three variables \code{t}, \code{x} and \code{y} for process \eqn{X_t}{X(t)} and \eqn{Y_t}{Y(t)}.}
  \item{corr}{the correlation structure of two Brownian motions W1(t) and W2(t); must be a real symmetric positive-definite square matrix of dimension 2.}
  \item{alpha, mu}{weight of the predictor-corrector scheme; the default \code{alpha = 0.5} and \code{mu = 0.5}.}
  \item{type}{if \code{type="ito"} simulation diffusion bridge of Itô type, else \code{type="str"} simulation diffusion bridge of Stratonovich type; the default \code{type="ito"}.}
  \item{method}{numerical methods of simulation, the default \code{method = "euler"}; see \code{\link{snssde2d}}.}
  \item{x, object}{an object inheriting from class \code{"bridgesde2d"}.}
  \item{at}{time between \code{t0} and \code{T}. Monte-Carlo statistics of the solution \eqn{(X_{t},Y_{t})}{(X(t),Y(t))} at time \code{at}. The default \code{at = T/2}.} 
  \item{digits}{integer, used for number formatting.} 
  \item{\dots}{potentially further arguments for (non-default) methods.}
}
\details{
The function \code{bridgesde2d} returns a \code{mts} of the diffusion bridge starting at \code{x} at time \code{t0} and ending at \code{y} at time \code{T}. W1(t) and W2(t) are two standard Brownian motion independent if \code{corr=NULL}.

The methods of approximation are classified according to their different properties. Mainly two criteria of optimality are used in the literature: the strong 
and the weak (orders of) convergence. The \code{method} of simulation can be one among: \code{Euler-Maruyama Order 0.5}, \code{Milstein Order 1}, \code{Milstein Second-Order}, 
\code{Predictor-Corrector method}, \code{Itô-Taylor Order 1.5}, \code{Heun Order 2} and \code{Runge-Kutta Order 1, 2 and 3}.

An overview of this package, see \code{\link[utils]{browseVignettes}('Sim.DiffProc')} for more informations.
}
\value{
\code{bridgesde2d} returns an object inheriting from \code{\link[base]{class}} \code{"bridgesde2d"}.
\item{X, Y }{an invisible \code{mts} (2-dim) object (X(t),Y(t)).}
\item{driftx, drifty }{drift coefficient of X(t) and Y(t).}
\item{diffx, diffy }{diffusion coefficient of X(t) and Y(t).}
\item{Cx, Cy}{indices of crossing realized of X(t) and Y(t).}
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
\code{\link{bridgesde1d}} for simulation of 1-dim SDE. 

\code{\link[sde]{DBridge}} in package "sde".
}
\examples{
## dX(t) = 4*(-1-X(t)) dt + 0.2 dW1(t)
## dY(t) = X(t) dt + 0 dW2(t)
## x01 = 0 , y01 = 0
## x02 = 0, y02 = 0 
## W1(t) and W2(t) two correlated Brownian motion with matrix Sigma=matrix(c(1,0.7,0.7,1),nrow=2)
set.seed(1234)

fx <- expression(4*(-1-x) , x)
gx <- expression(0.2 , 0)
Sigma= matrix(c(1,0.7,0.7,1),nrow=2)
res <- bridgesde2d(drift=fx,diffusion=gx,Dt=0.005,M=500,corr=Sigma)
res
summary(res) ## Monte-Carlo statistics at time T/2=2.5
summary(res,at=1) ## Monte-Carlo statistics at time 1
summary(res,at=4) ## Monte-Carlo statistics at time 4
##
plot(res,type="n")
lines(time(res),apply(res$X,1,mean),col=3,lwd=2)
lines(time(res),apply(res$Y,1,mean),col=4,lwd=2)
legend("topright",c(expression(E(X[t])),expression(E(Y[t]))),lty=1,inset = .7,col=c(3,4))
##
plot2d(res)
}
\keyword{sde}
\keyword{ts}
\keyword{mts}