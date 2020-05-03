

\name{bridgesde3d}
\alias{bridgesde3d}
\alias{bridgesde3d.default}
\alias{print.bridgesde3d}
\alias{time.bridgesde3d}
\alias{mean.bridgesde3d}
\alias{Median.bridgesde3d}
\alias{Mode.bridgesde3d}
\alias{quantile.bridgesde3d}
\alias{kurtosis.bridgesde3d}
\alias{skewness.bridgesde3d}
\alias{moment.bridgesde3d}
\alias{summary.bridgesde3d}
\alias{cv.bridgesde3d}
\alias{max.bridgesde3d}
\alias{min.bridgesde3d}
\alias{bconfint.bridgesde3d}
\alias{plot.bridgesde3d}
\alias{points.bridgesde3d}
\alias{lines.bridgesde3d}
\alias{plot3D.bridgesde3d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulation of 3-D Bridge SDE's}
\description{The (S3) generic function \code{bridgesde3d} for simulation of 3-dim bridge stochastic differential equations,Itô or Stratonovich type, with different methods.}
\usage{
bridgesde3d(N, \dots)
\method{bridgesde3d}{default}(N=1000,M=1, x0=c(0,0,0), 
   y=c(0,0,0), t0 = 0, T = 1, Dt, drift, diffusion, corr = NULL, 
   alpha = 0.5, mu = 0.5,type = c("ito", "str"), method = 
   c("euler", "milstein","predcorr","smilstein", "taylor", 
	"heun","rk1", "rk2", "rk3"), \dots)

\method{summary}{bridgesde3d}(object, at,
     digits=NULL, \dots)								  
\method{time}{bridgesde3d}(x, \dots)
\method{mean}{bridgesde3d}(x, at, \dots)
\method{Median}{bridgesde3d}(x, at, \dots)
\method{Mode}{bridgesde3d}(x, at, \dots)
\method{quantile}{bridgesde3d}(x, at, \dots)
\method{kurtosis}{bridgesde3d}(x, at, \dots)
\method{skewness}{bridgesde3d}(x, at, \dots)
\method{min}{bridgesde3d}(x, at, \dots)
\method{max}{bridgesde3d}(x, at, \dots)
\method{moment}{bridgesde3d}(x, at, \dots)
\method{cv}{bridgesde3d}(x, at, \dots)
\method{bconfint}{bridgesde3d}(x, at, \dots)

\method{plot}{bridgesde3d}(x, \dots)
\method{lines}{bridgesde3d}(x, \dots)
\method{points}{bridgesde3d}(x, \dots)	
\method{plot3D}{bridgesde3d}(x, display = c("persp","rgl"), \dots)						  
}
%- maybe also 'usage' for other xs documented here.
\arguments{
  \item{N}{number of simulation steps.}
  \item{M}{number of trajectories.}
  \item{x0}{initial value (numeric vector of length 3) of the process \eqn{X_t}{X(t)}, \eqn{Y_t}{Y(t)} and \eqn{Z_t}{Z(t)} at time \eqn{t_0}{t0}.}
  \item{y}{terminal value (numeric vector of length 3) of the process \eqn{X_t}{X(t)}, \eqn{Y_t}{Y(t)} and \eqn{Z_t}{Z(t)} at time \eqn{T}{T}.}
  \item{t0}{initial time.}
  \item{T}{final time.} 
  \item{Dt}{time step of the simulation (discretization). If it is \code{\link[base]{missing}} a default \eqn{\Delta t = \frac{T-t_{0}}{N}}{Dt = (T-t0)/N}.}
  \item{drift}{drift coefficient: an \code{\link[base]{expression}} of four variables \code{t}, \code{x}, \code{y} and \code{z} for process \eqn{X_t}{X(t)}, \eqn{Y_t}{Y(t)} and \eqn{Z_t}{Z(t)}.}
  \item{diffusion}{diffusion coefficient: an \code{\link[base]{expression}} of four variables \code{t}, \code{x}, \code{y} and \code{z} for process \eqn{X_t}{X(t)}, \eqn{Y_t}{Y(t)} and \eqn{Z_t}{Z(t)}.}
  \item{corr}{the correlation structure of three Brownian motions W1(t), W2(t) and W3(t); must be a real symmetric positive-definite square matrix of dimension 3.}
  \item{alpha}{weight \code{alpha} of the predictor-corrector scheme; the default \code{alpha = 0.5}.}
  \item{mu}{weight \code{mu} of the predictor-corrector scheme; the default \code{mu = 0.5}.}
  \item{type}{if \code{type="ito"} simulation diffusion bridge of Itô type, else \code{type="str"} simulation diffusion bridge of Stratonovich type; the default \code{type="ito"}.}
  \item{method}{numerical methods of simulation, the default \code{method = "euler"}; see \code{\link{snssde3d}}.}
  \item{x, object}{an object inheriting from class \code{"bridgesde3d"}.}
  \item{at}{time between \code{t0} and \code{T}. Monte-Carlo statistics of the solution \eqn{(X_{t},Y_{t},Z_{t})}{(X(t),Y(t),Z(t))} at time \code{at}. The default \code{at = T/2}.} 
  \item{digits}{integer, used for number formatting.} 
  \item{display}{\code{"persp"} perspective and \code{"rgl"} plots.} 
  \item{\dots}{potentially further arguments for (non-default) methods.}
}
\details{
The function \code{bridgesde3d} returns a \code{mts} of the diffusion bridge starting at \code{x} at time \code{t0} and ending at \code{y} at time \code{T}. W1(t), W2(t) and W3(t) three standard Brownian motion independent if \code{corr=NULL}.

The methods of approximation are classified according to their different properties. Mainly two criteria of optimality are used in the literature: the strong 
and the weak (orders of) convergence. The \code{method} of simulation can be one among: \code{Euler-Maruyama Order 0.5}, \code{Milstein Order 1}, \code{Milstein Second-Order}, 
\code{Predictor-Corrector method}, \code{Itô-Taylor Order 1.5}, \code{Heun Order 2} and \code{Runge-Kutta Order 1, 2 and 3}. 

An overview of this package, see \code{\link[utils]{browseVignettes}('Sim.DiffProc')} for more informations.
}
\value{
\code{bridgesde3d} returns an object inheriting from \code{\link[base]{class}} \code{"bridgesde3d"}.
\item{X, Y, Z }{an invisible \code{mts} (3-dim) object (X(t),Y(t),Z(t)).}
\item{driftx, drifty, driftz }{drift coefficient of X(t), Y(t) and Z(t).}
\item{diffx, diffy, diffz }{diffusion coefficient of X(t), Y(t) and Z(t).}
\item{Cx, Cy, Cz}{indices of crossing realized of X(t), Y(t)) and Z(t).}
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
\code{\link{bridgesde1d}} for simulation of 1-dim SDE. \code{\link[sde]{DBridge}} in package "sde".

\code{\link{bridgesde2d}} for simulation of 2-dim SDE.
}
\examples{
## dX(t) = 4*(-1-X(t))*Y(t) dt + 0.2 * dW1(t) ; x01 = 0 and y01 = 0
## dY(t) = 4*(1-Y(t)) *X(t) dt + 0.2 * dW2(t) ; x02 = -1 and y02 = -2
## dZ(t) = 4*(1-Z(t)) *Y(t) dt + 0.2 * dW3(t) ; x03 = 0.5 and y03 = 0.5       
## W1(t), W2(t) and W3(t) are three correlated Brownian motions with Sigma
set.seed(1234)

fx <- expression(4*(-1-x)*y, 4*(1-y)*x, 4*(1-z)*y)
gx <- rep(expression(0.2),3)
# correlation matrix
Sigma <-matrix(c(1,0.3,0.5,0.3,1,0.2,0.5,0.2,1),nrow=3,ncol=3) 

res <- bridgesde3d(x0=c(0,-1,0.5),y=c(0,-2,0.5),drift=fx,diffusion=gx,corr=Sigma,M=200)
res
summary(res) ## Monte-Carlo statistics at time T/2=0.5
summary(res,at=0.25) ## Monte-Carlo statistics at time 0.25
summary(res,at=0.75) ## Monte-Carlo statistics at time 0.75
##
plot(res,type="n")
lines(time(res),apply(res$X,1,mean),col=3,lwd=2)
lines(time(res),apply(res$Y,1,mean),col=4,lwd=2)
lines(time(res),apply(res$Z,1,mean),col=5,lwd=2)
legend("topleft",c(expression(E(X[t])),expression(E(Y[t])),
       expression(E(Z[t]))),lty=1,inset = .01,col=c(3,4,5))
##
plot3D(res,display = "persp",main="3-dim bridge sde")
}
\keyword{sde}
\keyword{ts}
\keyword{mts}
