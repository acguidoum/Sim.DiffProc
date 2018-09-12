

\name{HWV}
\alias{HWV}
\alias{OU}
\alias{HWV.default}
\alias{OU.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hull-White/Vasicek, Ornstein-Uhlenbeck process
}
\description{
The (S3) generic function for simulation of Hull-White/Vasicek or gaussian diffusion models, and Ornstein-Uhlenbeck process.
}
\usage{
HWV(N, \dots)
OU(N, \dots)

\method{HWV}{default}(N = 100, M = 1, x0 = 2, t0 = 0, T = 1, Dt = NULL, mu = 4, theta = 1,
   sigma = 0.1, \dots)
\method{OU}{default}(N =100,M=1,x0=2,t0=0,T=1,Dt = NULL,mu=4,sigma=0.2, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{number of simulation steps.}
  \item{M}{number of trajectories.}
  \item{x0}{initial value of the process at time \eqn{t_{0}}{\code{t0}}.}
  \item{t0}{initial time.}
  \item{T}{final time.} 
  \item{Dt}{time step of the simulation (discretization). If it is \code{\link[base]{missing}} a default \eqn{\Delta t = \frac{T-t_{0}}{N}}{Dt = (T-t0)/N}.}
  \item{mu}{parameter of the \code{HWV} and \code{OU}; see details.}
  \item{theta}{parameter of the \code{HWV}; see details.}
  \item{sigma}{the volatility of the \code{HWV} and \code{OU}.}
  \item{\dots}{potentially further arguments for (non-default) methods.}
}
\details{
The function \code{HWV} returns a trajectory of the \bold{Hull-White/Vasicek process} starting at \eqn{x_{0}}{x0} at time \eqn{t_{0}}{t0}; 
i.e., the diffusion process solution of stochastic differential equation: 
\deqn{dX_{t}= \mu ( \theta -X_{t}) dt + \sigma dW_{t}}{ dX(t) = mu *( theta- X(t)) dt + sigma dW(t)}

The function \code{OU} returns a trajectory of the \bold{Ornstein-Uhlenbeck} starting at \eqn{x_{0}}{x0} at time \eqn{t_{0}}{t0}; 
i.e., the diffusion process solution of stochastic differential equation: 
\deqn{dX_{t}= -\mu X_{t} dt + \sigma dW_{t}}{ dX(t) = -mu * X(t) dt + sigma dW(t)}

Constraints: \eqn{\mu , \sigma >0}{mu, sigma >0}.

Please note that the process is stationary only if \eqn{\mu >0}{mu >0}.
}
\value{
\item{X}{an visible \code{ts} object.}
}
\references{
   Vasicek, O. (1977). 
   An Equilibrium Characterization of the Term Structure.
   \emph{Journal of Financial Economics},  5, 177--188. 
}
\author{
A.C. Guidoum, K. Boukhetala.
}

\seealso{
\code{\link[sde]{rcOU}} and \code{\link[sde]{rsOU}} for conditional and stationary law of Vasicek process are available in "sde".
}
\examples{
## Hull-White/Vasicek Models
## dX(t) = 4 * (2.5 - X(t)) * dt + 1 *dW(t), X0=10
set.seed(1234)

X <- HWV(N=1000,M=10,mu = 4, theta = 2.5,sigma = 1,x0=10)
plot(X,plot.type="single")
lines(as.vector(time(X)),rowMeans(X),col="red")

## Ornstein-Uhlenbeck Process
## dX(t) = -4 * X(t) * dt + 1 *dW(t) , X0=2
set.seed(1234)

X <- OU(N=1000,M=10,mu = 4,sigma = 1,x0=10)
plot(X,plot.type="single")
lines(as.vector(time(X)),rowMeans(X),col="red")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{sde}
\keyword{ts}
