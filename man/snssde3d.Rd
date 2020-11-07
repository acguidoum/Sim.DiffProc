

\name{snssde3d}
\alias{snssde3d}
\alias{snssde3d.default}
\alias{summary.snssde3d}
\alias{print.snssde3d}
\alias{time.snssde3d}
\alias{mean.snssde3d}
\alias{Median.snssde3d}
\alias{Mode.snssde3d}
\alias{quantile.snssde3d}
\alias{kurtosis.snssde3d}
\alias{skewness.snssde3d}
\alias{moment.snssde3d}
\alias{cv.snssde3d}
\alias{max.snssde3d}
\alias{min.snssde3d}
\alias{bconfint.snssde3d}
\alias{plot.snssde3d}
\alias{points.snssde3d}
\alias{lines.snssde3d}
\alias{plot3D.snssde3d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulation of 3-D Stochastic Differential Equation}
\description{The (S3) generic function \code{snssde3d} of simulation of solutions to 3-dim stochastic differential equations of Itô or Stratonovich type, with different methods.
}
\usage{
snssde3d(N, \dots)
\method{snssde3d}{default}(N = 1000, M =1, x0=c(0,0,0), t0 = 0, T = 1, Dt, 
   drift, diffusion, corr = NULL, alpha = 0.5, mu = 0.5, type = c("ito", "str"), 
   method = c("euler", "milstein","predcorr", "smilstein", "taylor", 
   "heun", "rk1", "rk2", "rk3"), \dots)	
   
\method{summary}{snssde3d}(object, at, digits=NULL,\dots)
\method{time}{snssde3d}(x, \dots)
\method{mean}{snssde3d}(x, at, \dots)
\method{Median}{snssde3d}(x, at, \dots)
\method{Mode}{snssde3d}(x, at, \dots)
\method{quantile}{snssde3d}(x, at, \dots)
\method{kurtosis}{snssde3d}(x, at, \dots)
\method{skewness}{snssde3d}(x, at, \dots)
\method{min}{snssde3d}(x, at, \dots)
\method{max}{snssde3d}(x, at, \dots)
\method{moment}{snssde3d}(x, at, \dots)
\method{cv}{snssde3d}(x, at, \dots)
\method{bconfint}{snssde3d}(x, at, \dots)

\method{plot}{snssde3d}(x, \dots)
\method{lines}{snssde3d}(x, \dots)
\method{points}{snssde3d}(x, \dots)
\method{plot3D}{snssde3d}(x, display = c("persp","rgl"), \dots)
}
%- maybe also 'usage' for other xs documented here.
\arguments{
  \item{N}{number of simulation steps.}
  \item{M}{number of trajectories.}
  \item{x0}{initial value of the process \eqn{X_{t}}{X(t)}, \eqn{Y_{t}}{Y(t)} and \eqn{Z_{t}}{Z(t)} at time \code{t0}.}
  \item{t0}{initial time.}
  \item{T}{ending time.} 
  \item{Dt}{time step of the simulation (discretization). If it is \code{\link[base]{missing}} a default \eqn{\Delta t = \frac{T-t_{0}}{N}}{Dt = (T-t0)/N}.}
  \item{drift}{drift coefficient: an \code{\link[base]{expression}} of four variables \code{t}, \code{x}, \code{y} and \code{z} for process \eqn{X_t}{X(t)}, \eqn{Y_t}{Y(t)} and \eqn{Z_t}{Z(t)}.}
  \item{diffusion}{diffusion coefficient: an \code{\link[base]{expression}} of four variables \code{t}, \code{x}, \code{y} and \code{z} for process \eqn{X_t}{X(t)}, \eqn{Y_t}{Y(t)} and \eqn{Z_t}{Z(t)}.}
  \item{corr}{the correlation structure of three Brownian motions W1(t), W2(t) and W3(t); must be a real symmetric positive-definite square matrix of dimension 3.}
  \item{alpha, mu }{weight of the predictor-corrector scheme; the default \code{alpha = 0.5} and \code{mu = 0.5}.}
  \item{type}{if \code{type="ito"} simulation sde of Itô type, else \code{type="str"} simulation sde of Stratonovich type; the default \code{type="ito"}.}
  \item{method}{numerical methods of simulation, the default \code{method = "euler"}.}
  \item{x, object}{an object inheriting from class \code{"snssde3d"}.}
  \item{at}{time between \code{t0} and \code{T}. Monte-Carlo statistics of the solutions \eqn{(X_{t},Y_{t},Z_{t})}{(X(t),Y(t),Z(t)} at time \code{at}. The default \code{at = T}.}  
  \item{digits}{integer, used for number formatting.}  
  \item{display}{\code{"persp"} perspective or \code{"rgl"} plots.} 
  \item{\dots}{potentially further arguments for (non-default) methods.}
}
\details{
The function \code{snssde3d} returns a \code{mts} x of length \code{N+1}; i.e. solution of the 3-dim sde \eqn{(X_{t},Y_{t},Z_{t})}{(X(t),Y(t),Z(t))}
of Ito or Stratonovich types; If \code{Dt} is not specified, then the best discretization \eqn{\Delta t = \frac{T-t_{0}}{N}}{Dt = (T-t0)/N}. 

The 3-dim Ito stochastic differential equation is:
\deqn{dX(t) = a(t,X(t),Y(t),Z(t)) dt + b(t,X(t),Y(t),Z(t)) dW_{1}(t)}{dX(t) = a(t,X(t),Y(t),Z(t))*dt + b(t,X(t),Y(t),Z(t))*dW1(t)}
\deqn{dY(t) = a(t,X(t),Y(t),Z(t)) dt + b(t,X(t),Y(t),Z(t)) dW_{2}(t)}{dY(t) = a(t,X(t),Y(t),Z(t))*dt + b(t,X(t),Y(t),Z(t))*dW2(t)}
\deqn{dZ(t) = a(t,X(t),Y(t),Z(t)) dt + b(t,X(t),Y(t),Z(t)) dW_{3}(t)}{dZ(t) = a(t,X(t),Y(t),Z(t))*dt + b(t,X(t),Y(t),Z(t))*dW3(t)}
3-dim Stratonovich sde :
\deqn{dX(t) = a(t,X(t),Y(t),Z(t)) dt + b(t,X(t),Y(t),Z(t)) \circ dW_{1}(t)}{dX(t) = a(t,X(t),Y(t),Z(t))*dt + b(t,X(t),Y(t),Z(t)) o dW1(t)}
\deqn{dY(t) = a(t,X(t),Y(t),Z(t)) dt + b(t,X(t),Y(t),Z(t)) \circ dW_{2}(t)}{dY(t) = a(t,X(t),Y(t),Z(t))*dt + b(t,X(t),Y(t),Z(t)) o dW2(t)}
\deqn{dZ(t) = a(t,X(t),Y(t),Z(t)) dt + b(t,X(t),Y(t),Z(t)) \circ dW_{3}(t)}{dZ(t) = a(t,X(t),Y(t),Z(t))*dt + b(t,X(t),Y(t),Z(t)) o dW3(t)}
\eqn{W_{1}(t), W_{2}(t), W_{3}(t)}{W1(t), W2(t) and W3(t)} three standard Brownian motion independent if \code{corr=NULL}.

In the correlation case, currently we can use only the Euler-Maruyama and Milstein scheme.

The methods of approximation are classified according to their different properties. Mainly two criteria of optimality are used in the literature: the strong 
and the weak (orders of) convergence. The \code{method} of simulation can be one among: \code{Euler-Maruyama Order 0.5}, \code{Milstein Order 1}, \code{Milstein Second-Order}, 
\code{Predictor-Corrector method}, \code{Itô-Taylor Order 1.5}, \code{Heun Order 2} and \code{Runge-Kutta Order 1, 2 and 3}. 

An overview of this package, see \code{\link[utils]{browseVignettes}('Sim.DiffProc')} for more informations.
}
\value{
\code{snssde3d} returns an object inheriting from \code{\link[base]{class}} \code{"snssde3d"}.
\item{X, Y, Z }{an invisible \code{mts} (3-dim) object (X(t),Y(t),Z(t)).}
\item{driftx, drifty, driftz }{drift coefficient of X(t), Y(t) and Z(t).}
\item{diffx, diffy, diffz }{diffusion coefficient of X(t), Y(t) and Z(t).}
\item{type }{type of sde.}
\item{method }{the numerical method used.}
}
\references{

   Guidoum AC, Boukhetala K (2020). 
   "Performing Parallel Monte Carlo and Moment Equations Methods for Itô and Stratonovich Stochastic Differential Systems: R Package Sim.DiffProc". 
   \emph{Journal of Statistical Software}, \bold{96}(2), 1--82. doi:10.18637/jss.v096.i02
   
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
\code{\link{snssde1d}} and \code{\link{snssde2d}} for 1- and 2-dim sde.

\code{sde.sim} in package "sde". \code{simulate} in package "yuima". 
}
\examples{

## Example 1: Ito sde
## dX(t) = (2*(Y(t)>0)-2*(Z(t)<=0)) dt + 0.2 * dW1(t) 
## dY(t) = -2*Y(t) dt + 0.2 * dW2(t) 
## dZ(t) = -2*Z(t) dt + 0.2 * dW3(t)        
## W1(t), W2(t) and W3(t) three independent Brownian motion
set.seed(1234)

fx <- expression(2*(y>0)-2*(z<=0) , -2*y, -2*z)
gx <- rep(expression(0.2),3)

mod3d1 <- snssde3d(x0=c(0,2,-2),drift=fx,diffusion=gx,M=500,Dt=0.003)
mod3d1
summary(mod3d1)
##
dev.new()
plot(mod3d1,type="n")
mx <- apply(mod3d1$X,1,mean)
my <- apply(mod3d1$Y,1,mean)
mz <- apply(mod3d1$Z,1,mean)
lines(time(mod3d1),mx,col=1)
lines(time(mod3d1),my,col=2)
lines(time(mod3d1),mz,col=3)
legend("topright",c(expression(E(X[t])),expression(E(Y[t])),
 expression(E(Z[t]))),lty=1,inset = .01,col=c(1,2,3),cex=0.95)
##
dev.new()
plot3D(mod3d1,display="persp") ## in space (O,X,Y,Z)

## Example 2: Stratonovich sde 
## dX(t) = Y(t)* dt  + 0.2 o dW3(t)         
## dY(t) = (4*( 1-X(t)^2 )* Y(t) - X(t))* dt + 0.2 o dW2(t)
## dZ(t) = (4*( 1-X(t)^2 )* Z(t) - X(t))* dt + 0.2 o dW3(t)
## W1(t), W2(t) and W3(t) are three correlated Brownian motions with Sigma

fx <- expression( y , (4*( 1-x^2 )* y - x), (4*( 1-x^2 )* z - x))
gx <- expression( 0.2 , 0.2, 0.2)
# correlation matrix
Sigma <-matrix(c(1,0.3,0.5,0.3,1,0.2,0.5,0.2,1),nrow=3,ncol=3) 

mod3d2 <- snssde3d(drift=fx,diffusion=gx,N=10000,T=100,type="str",corr=Sigma)
mod3d2
##
dev.new()
plot(mod3d2,pos=2)
##
dev.new()
plot(mod3d2,union = FALSE)
##
dev.new()
plot3D(mod3d2,display="persp") ## in space (O,X,Y,Z)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{sde}
\keyword{ts}
\keyword{mts}
