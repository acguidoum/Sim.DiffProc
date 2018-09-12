

\name{BM}
\alias{ABM}
\alias{BB}
\alias{BM}
\alias{GBM}
\alias{ABM.default}
\alias{BB.default}
\alias{BM.default}
\alias{GBM.default}
\title{
Brownian motion, Brownian bridge, geometric Brownian motion, and arithmetic Brownian motion simulators
}
\description{
The (S3) generic function for simulation of brownian motion, brownian bridge, geometric brownian motion, and arithmetic brownian motion.
}
\usage{
BM(N, \dots)
BB(N, \dots)
GBM(N, \dots)
ABM(N, \dots)

\method{BM}{default}(N =1000,M=1,x0=0,t0=0,T=1,Dt=NULL, \dots)
\method{BB}{default}(N =1000,M=1,x0=0,y=0,t0=0,T=1,Dt=NULL, \dots)
\method{GBM}{default}(N =1000,M=1,x0=1,t0=0,T=1,Dt=NULL,theta=1,sigma=1, \dots)
\method{ABM}{default}(N =1000,M=1,x0=0,t0=0,T=1,Dt=NULL,theta=1,sigma=1, \dots)
}
\arguments{
  \item{N}{number of simulation steps.}
  \item{M}{number of trajectories.}
  \item{x0}{initial value of the process at time \eqn{t_{0}}{\code{t0}}.}
  \item{y}{terminal value of the process at time \eqn{T}{\code{T}} of the \code{BB}.}
  \item{t0}{initial time.}
  \item{T}{final time.} 
  \item{Dt}{time step of the simulation (discretization). If it is \code{\link[base]{NULL}} a default \eqn{\Delta t = \frac{T-t_{0}}{N}}{Dt = (T-t0)/N}.}
  \item{theta}{the interest rate of the \code{ABM} and \code{GBM}.}
  \item{sigma}{the volatility of the \code{ABM} and \code{GBM}.}
  \item{\dots}{potentially further arguments for (non-default) methods.}
}
\details{
The function \code{BM} returns a trajectory of the \bold{standard Brownian motion} (Wiener process) in the time interval \eqn{[t_{0},T]}{[t0,T]}. Indeed, for \eqn{W(dt)}{W(dt)} it holds true that 
\eqn{W(dt) \rightarrow W(dt) - W(0) \rightarrow \mathcal{N}(0,dt)}{W(dt) = W(dt) - W(0) -> N(0,dt) -> sqrt(dt) * N(0,1)}, where \eqn{\mathcal{N}(0,1)}{N(0,1)} is normal distribution 
\code{\link[stats]{Normal}}.

The function \code{BB} returns a trajectory of the \bold{Brownian bridge} starting at \eqn{x_{0}}{x0} at time \eqn{t_{0}}{t0} and ending
at \eqn{y}{y} at time \eqn{T}{T}; i.e., the diffusion process solution of stochastic differential equation: 
\deqn{dX_{t}= \frac{y-X_{t}}{T-t} dt + dW_{t}}{ dX(t) = ((y-X(t))/(T-t)) dt + dW(t)}

The function \code{GBM} returns a trajectory of the \bold{geometric Brownian motion} starting at \eqn{x_{0}}{x0} at time \eqn{t_{0}}{t0}; 
i.e., the diffusion process solution of stochastic differential equation: 
\deqn{dX_{t}= \theta X_{t} dt + \sigma X_{t} dW_{t}}{ dX(t) = theta X(t) dt + sigma X(t) dW(t)}

The function \code{ABM} returns a trajectory of the \bold{arithmetic Brownian motion} starting at \eqn{x_{0}}{x0} at time \eqn{t_{0}}{t0}; 
i.e.,; the diffusion process solution of stochastic differential equation: 
\deqn{dX_{t}= \theta dt + \sigma dW_{t}}{ dX(t) = theta dt + sigma dW(t)}
}
\value{
\item{X}{an visible \code{ts} object.}
}
\references{   
   Allen, E. (2007).
   \emph{Modeling with Ito stochastic differential equations.}
   Springer-Verlag, New York.
   
   Jedrzejewski, F. (2009).
   \emph{Modeles aleatoires et physique probabiliste.}
   Springer-Verlag, New York.

   Henderson, D and Plaschko, P. (2006). 
   \emph{Stochastic differential equations in science and engineering.}
   World Scientific.   
}
\author{
A.C. Guidoum, K. Boukhetala.
}

\seealso{
This functions \code{\link[sde]{BM}}, \code{\link[sde]{BBridge}} and \code{\link[sde]{GBM}} are available in other packages such as "sde".
}
\examples{

op <- par(mfrow = c(2, 2))

## Brownian motion
set.seed(1234)
X <- BM(M = 100)
plot(X,plot.type="single")
lines(as.vector(time(X)),rowMeans(X),col="red")

## Brownian bridge
set.seed(1234)
X <- BB(M =100)
plot(X,plot.type="single")
lines(as.vector(time(X)),rowMeans(X),col="red")

## Geometric Brownian motion
set.seed(1234)
X <- GBM(M = 100)
plot(X,plot.type="single")
lines(as.vector(time(X)),rowMeans(X),col="red")

## Arithmetic Brownian motion
set.seed(1234)
X <- ABM(M = 100)
plot(X,plot.type="single")
lines(as.vector(time(X)),rowMeans(X),col="red")

par(op)
}
\keyword{BM}
\keyword{sde}
\keyword{ts}
