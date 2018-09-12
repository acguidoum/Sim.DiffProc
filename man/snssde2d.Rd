

\name{snssde2d}
\alias{snssde2d}
\alias{snssde2d.default}
\alias{summary.snssde2d}
\alias{print.snssde2d}
\alias{time.snssde2d}
\alias{mean.snssde2d}
\alias{Median.snssde2d}
\alias{Mode.snssde2d}
\alias{quantile.snssde2d}
\alias{kurtosis.snssde2d}
\alias{skewness.snssde2d}
\alias{moment.snssde2d}
\alias{cv.snssde2d}
\alias{max.snssde2d}
\alias{min.snssde2d}
\alias{bconfint.snssde2d}
\alias{plot.snssde2d}
\alias{points.snssde2d}
\alias{lines.snssde2d}
\alias{plot2d.snssde2d}
\alias{points2d.snssde2d}
\alias{lines2d.snssde2d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulation of 2-D Stochastic Differential Equation}
\description{The (S3) generic function \code{snssde2d} of simulation of solutions to 2-dim stochastic differential equations of Itô or Stratonovich type, with different methods.
}
\usage{
snssde2d(N, \dots)
\method{snssde2d}{default}(N = 1000, M = 1, x0 = c(0,0),t0 = 0, T = 1, Dt=NULL, 
   drift,diffusion, alpha = 0.5, mu = 0.5, type = c("ito", "str"), 
   method = c("euler", "milstein", "predcorr", "smilstein", "taylor", 
   "heun", "rk1", "rk2","rk3"), \dots)	
   
\method{summary}{snssde2d}(object, at, digits=NULL,\dots)
\method{time}{snssde2d}(x, \dots)
\method{mean}{snssde2d}(x, at, \dots)
\method{Median}{snssde2d}(x, at, \dots)
\method{Mode}{snssde2d}(x, at, \dots)
\method{quantile}{snssde2d}(x, at, \dots)
\method{kurtosis}{snssde2d}(x, at, \dots)
\method{skewness}{snssde2d}(x, at, \dots)
\method{min}{snssde2d}(x, at, \dots)
\method{max}{snssde2d}(x, at, \dots)
\method{moment}{snssde2d}(x, at, \dots)
\method{cv}{snssde2d}(x, at, \dots)
\method{bconfint}{snssde2d}(x, at, \dots)

\method{plot}{snssde2d}(x, \dots)
\method{lines}{snssde2d}(x, \dots)
\method{points}{snssde2d}(x, \dots)
\method{plot2d}{snssde2d}(x, \dots)
\method{lines2d}{snssde2d}(x, \dots)
\method{points2d}{snssde2d}(x, \dots)
}
%- maybe also 'usage' for other xs documented here.
\arguments{
  \item{N}{number of simulation steps.}
  \item{M}{number of trajectories (Monte-Carlo).}
  \item{x0}{initial values x0=(x,y) of the process \eqn{X_{t}}{X(t)} and \eqn{Y_{t}}{Y(t)} at time \code{t0}.}
  \item{t0}{initial time.}
  \item{T}{ending time.} 
  \item{Dt}{time step of the simulation (discretization). If it is \code{\link[base]{NULL}} a default \eqn{\Delta t = \frac{T-t_{0}}{N}}{Dt = (T-t0)/N}.}
  \item{drift}{drift coefficient: an \code{\link[base]{expression}} of three variables \code{t}, \code{x} and \code{y} for process \eqn{X_t}{X(t)} and \eqn{Y_t}{Y(t)}.}
  \item{diffusion}{diffusion coefficient: an \code{\link[base]{expression}} of three variables \code{t}, \code{x} and \code{y} for process \eqn{X_t}{X(t)} and \eqn{Y_t}{Y(t)}.}
  \item{alpha, mu}{weight of the predictor-corrector scheme; the default \code{alpha = 0.5} and \code{mu = 0.5}.}
  \item{type}{if \code{type="ito"} simulation sde of Itô type, else \code{type="str"} simulation sde of Stratonovich type; the default \code{type="ito"}.}
  \item{method}{numerical methods of simulation, the default \code{method = "euler"}.}
  \item{x, object}{an object inheriting from class \code{"snssde2d"}.} 
  \item{at}{time between \code{t0} and \code{T}. Monte-Carlo statistics of the solutions \eqn{(X_{t},Y_{t})}{(X(t),Y(t))} at time \code{at}. The default \code{at = T}.}    
  \item{digits}{integer, used for number formatting.}  
  \item{\dots}{potentially further arguments for (non-default) methods.}
}
\details{
The function \code{snssde2d} returns a \code{mts} x of length \code{N+1}; i.e. solution of the 2-dim sde \eqn{(X_{t},Y_{t})}{(X(t),Y(t))}
of Ito or Stratonovich types; If \code{Dt} is not specified, then the best discretization \eqn{\Delta t = \frac{T-t_{0}}{N}}{Dt = (T-t0)/N}. 

The 2-dim Ito stochastic differential equation is: 
\deqn{dX(t) = a(t,X(t),Y(t)) dt + b(t,X(t),Y(t)) dW_{1}(t)}{dX(t) = a(t,X(t),Y(t))*dt + b(t,X(t),Y(t))*dW1(t)}
\deqn{dY(t) = a(t,X(t),Y(t)) dt + b(t,X(t),Y(t)) dW_{2}(t)}{dY(t) = a(t,X(t),Y(t))*dt + b(t,X(t),Y(t))*dW2(t)}
2-dim Stratonovich sde :
\deqn{dX(t) = a(t,X(t),Y(t)) dt + b(t,X(t),Y(t)) \circ dW_{1}(t)}{dX(t) = a(t,X(t),Y(t))*dt + b(t,X(t),Y(t)) o dW1(t)}
\deqn{dY(t) = a(t,X(t),Y(t)) dt + b(t,X(t),Y(t)) \circ dW_{2}(t)}{dY(t) = a(t,X(t),Y(t))*dt + b(t,X(t),Y(t)) o dW2(t)}
\eqn{W_{1}(t), W_{2}(t)}{W1(t) and W2(t)} two standard Brownian motion independent.

The methods of approximation are classified according to their different properties. Mainly two criteria of optimality are used in the literature: the strong 
and the weak (orders of) convergence. The \code{method} of simulation can be one among: \code{Euler-Maruyama Order 0.5}, \code{Milstein Order 1}, \code{Milstein Second-Order}, 
\code{Predictor-Corrector method}, \code{Itô-Taylor Order 1.5}, \code{Heun Order 2} and \code{Runge-Kutta Order 1, 2 and 3}. 

An overview of this package, see \code{\link[utils]{browseVignettes}('Sim.DiffProc')} for more informations.
}
\value{
\code{snssde2d} returns an object inheriting from \code{\link[base]{class}} \code{"snssde2d"}.
\item{X, Y }{an invisible \code{mts} (2-dim) object (X(t),Y(t)).}
\item{driftx, drifty }{drift coefficient of X(t) and Y(t).}
\item{diffx, diffy }{diffusion coefficient of X(t) and Y(t).}
\item{type }{type of sde.}
\item{method }{the numerical method used.}
}
\references{
  Friedman, A. (1975). 
  \emph{Stochastic differential equations and applications}.
  Volume 1, ACADEMIC PRESS.
  
  Henderson, D. and Plaschko,P. (2006). 
  \emph{Stochastic differential equations in science and engineering}.
  World Scientific.
  
  Allen, E. (2007).
  \emph{Modeling with Ito stochastic differential equations}. 
  Springer-Verlag.

  Jedrzejewski, F. (2009). 
  \emph{Modeles aleatoires et physique probabiliste}. 
  Springer-Verlag.
  
  Iacus, S.M. (2008). 
  \emph{Simulation and inference for stochastic differential equations: with R examples}. 
  Springer-Verlag, New York.
    
  Kloeden, P.E, and Platen, E. (1989). 
  A survey of numerical methods for stochastic differential equations. 
  \emph{Stochastic Hydrology and Hydraulics}, \bold{3}, 155--178.

  Kloeden, P.E, and Platen, E. (1991a).
  Relations between multiple ito and stratonovich integrals.
  \emph{Stochastic Analysis and Applications}, \bold{9}(3), 311--321.

  Kloeden, P.E, and Platen, E. (1991b).
  Stratonovich and ito stochastic taylor expansions. 
  \emph{Mathematische Nachrichten}, \bold{151}, 33--50.

  Kloeden, P.E, and Platen, E. (1995).
  \emph{Numerical Solution of Stochastic Differential Equations}. 
  Springer-Verlag, New York.

  Oksendal, B. (2000).
  \emph{Stochastic Differential Equations: An Introduction with Applications}.
  5th edn. Springer-Verlag, Berlin.

  Platen, E. (1980).
  Weak convergence of approximations of ito integral equations. 
  \emph{Z Angew Math Mech}. \bold{60}, 609--614.

  Platen, E. and Bruti-Liberati, N. (2010).
  \emph{Numerical Solution of Stochastic Differential Equations with Jumps in Finance.} 
  Springer-Verlag, New York

  Saito, Y, and Mitsui, T. (1993).
  Simulation of Stochastic Differential Equations. 
  \emph{The Annals of the Institute of Statistical Mathematics}, \bold{3}, 419--432.
}
\author{
A.C. Guidoum, K. Boukhetala.
}
\seealso{
\code{\link{snssde3d}} for 3-dim sde.

\code{simulate} in package "yuima". 
}
\examples{

## Example 1: Ito sde
## dX(t) = 4*(-1-X(t))*Y(t) dt + 0.2 dW1(t)
## dY(t) = 4*(1-Y(t))*X(t) dt + 0.2 dW2(t)
set.seed(1234)

fx <- expression(4*(-1-x)*y , 4*(1-y)*x )
gx <- expression(0.25*y,0.2*x)

mod2d1 <- snssde2d(drift=fx,diffusion=gx,x0=c(x0=1,y0=-1),M=1000)
mod2d1
summary(mod2d1)
##
dev.new()
plot(mod2d1,type="n")
mx <- apply(mod2d1$X,1,mean)
my <- apply(mod2d1$Y,1,mean)
lines(time(mod2d1),mx,col=1)
lines(time(mod2d1),my,col=2)
legend("topright",c(expression(E(X[t])),expression(E(Y[t]))),lty=1,inset = .01,col=c(1,2),cex=0.95)
##
dev.new()
plot2d(mod2d1) ## in plane (O,X,Y)
lines(my~mx,col=2)



## Example 2: Stratonovich sde
## dX(t) = Y(t) dt + 0 o dW1(t)
## dY(t) = (4*(1-X(t)^2)*Y(t) - X(t) ) dt + 0.2 o dW2(t)
set.seed(1234)

fx <- expression( y , (4*( 1-x^2 )* y - x))
gx <- expression( 0 , 0.2)

mod2d2 <- snssde2d(drift=fx,diffusion=gx,type="str",T=100,N=10000)
mod2d2
plot(mod2d2,pos=2)
dev.new()
plot(mod2d2,union = FALSE)
dev.new()
plot2d(mod2d2,type="n") ## in plane (O,X,Y)
points2d(mod2d2,col=rgb(0,100,0,50,maxColorValue=255), pch=16)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{sde}
\keyword{ts}
\keyword{mts}
