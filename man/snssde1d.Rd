

\name{snssde1d}
\alias{snssde1d}
\alias{snssde1d.default}
\alias{summary.snssde1d}
\alias{print.snssde1d}
\alias{time.snssde1d}
\alias{mean.snssde1d}
\alias{Median.snssde1d}
\alias{Mode.snssde1d}
\alias{quantile.snssde1d}
\alias{kurtosis.snssde1d}
\alias{skewness.snssde1d}
\alias{moment.snssde1d}
\alias{cv.snssde1d}
\alias{max.snssde1d}
\alias{min.snssde1d}
\alias{bconfint.snssde1d}
\alias{plot.snssde1d}
\alias{points.snssde1d}
\alias{lines.snssde1d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulation of 1-D Stochastic Differential Equation}
\description{The (S3) generic function \code{snssde1d} of simulation of solution to 1-dim stochastic differential equation of Itô or Stratonovich type, with different methods.}
\usage{
snssde1d(N, \dots)
\method{snssde1d}{default}(N = 1000, M = 1, x0 = 0, t0 = 0, T = 1, Dt, 
   drift, diffusion, alpha = 0.5, mu = 0.5, type = c("ito", "str"), 
   method = c("euler", "milstein", "predcorr", "smilstein", "taylor", 
   "heun", "rk1", "rk2", "rk3"), \dots)	
   
\method{summary}{snssde1d}(object, at ,digits=NULL, \dots)
\method{time}{snssde1d}(x, \dots)
\method{mean}{snssde1d}(x, at, \dots)
\method{Median}{snssde1d}(x, at, \dots)
\method{Mode}{snssde1d}(x, at, \dots)
\method{quantile}{snssde1d}(x, at, \dots)
\method{kurtosis}{snssde1d}(x, at, \dots)
\method{min}{snssde1d}(x, at, \dots)
\method{max}{snssde1d}(x, at, \dots)
\method{skewness}{snssde1d}(x, at, \dots)
\method{moment}{snssde1d}(x, at,  \dots)
\method{cv}{snssde1d}(x, at,  \dots)
\method{bconfint}{snssde1d}(x, at,  \dots)

\method{plot}{snssde1d}(x, \dots)
\method{lines}{snssde1d}(x, \dots)
\method{points}{snssde1d}(x, \dots)		
}
%- maybe also 'usage' for other xs documented here.
\arguments{
  \item{N}{number of simulation steps.}
  \item{M}{number of trajectories (Monte-Carlo).}
  \item{x0}{initial value of the process at time \code{t0}.}
  \item{t0}{initial time.}
  \item{T}{ending time.} 
  \item{Dt}{time step of the simulation (discretization). If it is \code{\link[base]{missing}} a default \eqn{\Delta t = \frac{T-t_{0}}{N}}{Dt = (T-t0)/N}.}
  \item{drift}{drift coefficient: an \code{\link[base]{expression}} of two variables \code{t} and \code{x}.}
  \item{diffusion}{diffusion coefficient: an \code{\link[base]{expression}} of two variables \code{t} and \code{x}.}
  \item{alpha, mu}{weight of the predictor-corrector scheme; the default \code{alpha = 0.5} and \code{mu = 0.5}.}
  \item{type}{if \code{type="ito"} simulation sde of Itô type, else \code{type="str"} simulation sde of Stratonovich type; the default \code{type="ito"}.}
  \item{method}{numerical methods of simulation, the default \code{method = "euler"}.}
  \item{x, object}{an object inheriting from class \code{"snssde1d"}.}
  \item{at}{time between \code{t0} and \code{T}. Monte-Carlo statistics of the solution \eqn{X_{t}}{X(t)} at time \code{at}. The default \code{at = T}.}   
  \item{digits}{integer, used for number formatting.}  
  \item{\dots}{potentially further arguments for (non-default) methods.}
}
\details{
The function \code{snssde1d} returns a \code{\link[stats]{ts}} x of length \code{N+1}; i.e. solution of the sde of Ito or 
Stratonovich types; If \code{Dt} is not specified, then the best discretization \eqn{\Delta t = \frac{T-t_{0}}{N}}{Dt = (T-t0)/N}. 

The Ito stochastic differential equation is: 
\deqn{dX(t) = a(t,X(t)) dt + b(t,X(t)) dW(t)}{dX(t) = a(t,X(t))*dt + b(t,X(t))*dW(t)}
Stratonovich sde :
\deqn{dX(t) = a(t,X(t)) dt + b(t,X(t)) \circ dW(t)}{dX(t) = a(t,X(t))*dt + b(t,X(t)) o dW(t)}

The methods of approximation are classified according to their different properties. Mainly two criteria of optimality are used in the literature: the strong 
and the weak (orders of) convergence. The \code{method} of simulation can be one among: \code{Euler-Maruyama Order 0.5}, \code{Milstein Order 1}, \code{Milstein Second-Order}, 
\code{Predictor-Corrector method}, \code{Itô-Taylor Order 1.5}, \code{Heun Order 2} and \code{Runge-Kutta Order 1, 2 and 3}. 

An overview of this package, see \code{\link[utils]{browseVignettes}('Sim.DiffProc')} for more informations.
}
\value{
\code{snssde1d} returns an object inheriting from \code{\link[base]{class}} \code{"snssde1d"}.
\item{X }{an invisible \code{\link[stats]{ts}} object.}
\item{drift }{drift coefficient.}
\item{diffusion }{diffusion coefficient.}
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
\code{\link{snssde2d}} and \code{\link{snssde3d}} for 2 and 3-dim sde.

\code{\link[sde]{sde.sim}} in package "sde". 

\code{simulate} in package "yuima".
}
\examples{

## Example 1: Ito sde
## dX(t) = 2*(3-X(t)) dt + 2*X(t) dW(t)
set.seed(1234)

f <- expression(2*(3-x) )
g <- expression(1)
mod1 <- snssde1d(drift=f,diffusion=g,M=4000,x0=10,Dt=0.01)
mod1
summary(mod1)
\dontrun{
plot(mod1)
lines(time(mod1),apply(mod1$X,1,mean),col=2,lwd=2)
lines(time(mod1),apply(mod1$X,1,bconfint,level=0.95)[1,],col=4,lwd=2)
lines(time(mod1),apply(mod1$X,1,bconfint,level=0.95)[2,],col=4,lwd=2)
legend("topright",c("mean path",paste("bound of", 95," percent confidence")),
       inset = .01,col=c(2,4),lwd=2,cex=0.8)
}
## Example 2: Stratonovich sde
## dX(t) = ((2-X(t))/(2-t)) dt + X(t) o dW(t)
set.seed(1234)

f <- expression((2-x)/(2-t))
g <- expression(x)
mod2 <- snssde1d(type="str",drift=f,diffusion=g,M=4000,x0=1, method="milstein")
mod2
summary(mod2,at = 0.25)
summary(mod2,at = 1)
\dontrun{
plot(mod2)
lines(time(mod2),apply(mod2$X,1,mean),col=2,lwd=2)
lines(time(mod2),apply(mod2$X,1,bconfint,level=0.95)[1,],col=4,lwd=2)
lines(time(mod2),apply(mod2$X,1,bconfint,level=0.95)[2,],col=4,lwd=2)
legend("topleft",c("mean path",paste("bound of", 95," percent confidence")),
       inset = .01,col=c(2,4),lwd=2,cex=0.8)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{sde}
\keyword{ts}
\keyword{mts}
