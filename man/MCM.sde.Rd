\name{MCM.sde}
\alias{MCM.sde}
\alias{MCM.sde.default}
\alias{print.MCM.sde}
\alias{plot.MCM.sde}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parallel Monte-Carlo Methods for SDE's
}
\description{
Generate \code{R} Monte-Carlo (version parallel) replicates of a statistic applied to SDE's (1,2 and 3 dim) for the two cases Ito and Stratonovich interpretations.
}
\usage{
MCM.sde(model, \dots)

\method{MCM.sde}{default}(model, statistic, R = 100, time, exact = NULL, 
         names = NULL, level = 0.95, parallel = c("no", "multicore", "snow"),
         ncpus = getOption("ncpus", 1L), cl = NULL, ...)

\method{plot}{MCM.sde}(x,index = 1,type=c("all","hist","qqplot","boxplot","CI"), \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
an object from class \code{\link{snssde1d}}, \code{\link{snssde2d}} and \code{\link{snssde3d}}.
}
  \item{statistic}{
a function which when applied to model returns a vector containing the statistic(s) of interest. 
}
  \item{R}{
the number of Monte-Carlo replicates. Usually this will be a single positive integer "R > 1".
}
  \item{time}{
the time when estimating the statistic(s) of interesttime between \code{t0} and \code{T}. The default \code{time = T}.
}
  \item{exact}{
a named list giving the exact statistic(s) if it exists otherwise \code{exact = NULL}. 
}
  \item{names}{
named the statistic(s) of interest. The default \code{names=c("mu1","mu2",...)}.
}
  \item{level}{
the confidence level(s) of the required interval(s).
}
  \item{parallel}{
the type of parallel operation to be used (if any). The default \code{parallel = "no"}.
}
  \item{ncpus}{
integer: number of processes to be used in parallel operation: typically one would chose this to the number of available CPUs.
}
  \item{cl}{
an optional parallel or snow cluster for use if \code{parallel = "snow"}.
}
  \item{x}{
an object inheriting from class \code{"MCM.sde"}.
}
  \item{index}{
the index of the variable of interest within the output of \code{"MCM.sde"}.
}
  \item{type}{
the type of plot of the Monte-Carlo estimation of the variable of interest. The default \code{type = "all"}.   
}
  \item{\dots}{
potentially further arguments for (non-default) methods.
}
}
\details{
We have here developed Monte-Carlo methods whose essence is the use of repeated experiments to evaluate a 
statistic(s) of interest in SDE's. For example estimation of moments as: mean, variance, covariance (and other as median, 
mode, quantile,...). With the standard error and the confidence interval for these estimators.

An overview of this package, see \code{\link[utils]{browseVignettes}('Sim.DiffProc')} for more informations.
}
\value{
The returned value is an object of class \code{"MCM.sde"}, containing the following components:

\item{mod}{The SDE's used (class: \code{\link{snssde1d}}, \code{\link{snssde2d}} and \code{\link{snssde3d}}).}
\item{dim}{Dimension of the model.}
\item{call}{The original call to \code{"MCM.sde"}.}
\item{Fn}{The function statistic as passed to \code{"MCM.sde"}.}
\item{ech}{A matrix with \code{sum(R)} column each of which is a Monte-Carlo replicate of the result of calling statistic.}
\item{time}{The time when estimating the statistic(s) of interest.}
\item{name}{named of statistic(s) of interest.}
\item{MC}{Table contains simulation results of statistic(s) of interest: Estimate, Bias (if exact available), Std.Error and Confidence interval.}
}
\references{

   Guidoum AC, Boukhetala K (2020). 
   "Performing Parallel Monte Carlo and Moment Equations Methods for Itô and Stratonovich Stochastic Differential Systems: R Package Sim.DiffProc". 
   \emph{Journal of Statistical Software}, \bold{96}(2), 1--82. doi:10.18637/jss.v096.i02
   
   Paul Glasserman (2003).
   \emph{Monte Carlo Methods in Financial Engineering}. 
   Springer-Verlag New York.

   Jun S. Liu (2004).
   \emph{Monte Carlo Strategies in Scientific Computing}. 
   Springer-Verlag New York.
   
   Christian Robert and George Casella (2010).
   \emph{Introducing Monte Carlo Methods with R}. 
   Springer-Verlag New York.
   
   Nick T. Thomopoulos (2013).
   \emph{Essentials of Monte Carlo Simulation: Statistical Methods for Building Simulation Models}. 
   Springer-Verlag New York.
   
   Q. Ethan McCallum and Stephen Weston (2011).
   \emph{Parallel R}. 
   O'Reilly Media, Inc.
}
\author{
A.C. Guidoum, K. Boukhetala.
}
\note{
When \code{parallel = "multicore"} is used are not available on Windows, \code{parallel = "snow"} is primarily intended to be used on multi-core Windows machine where \code{parallel = "multicore"} is not available. For more details see Q.E.McCallum and S.Weston (2011).
}

\seealso{
\code{\link{MEM.sde}} moment equations methods for SDE's.
}
\examples{
## Example 1 : (1 dim)
## dX(t) = 3*(1-X(t)) dt + 0.5 * dW(t), X(0)=5, t in [0,10]
## set the model 1d
f <- expression(3*(1-x));g <- expression(0.5)
mod1d <-  snssde1d(drift=f,diffusion=g,x0=5,T=10,M=50)

## function of the statistic(s) of interest.
sde.fun1d <- function(data, i){
  d <- data[i, ]
  return(c(mean(d),Mode(d),var(d)))
}

mc.sde1d = MCM.sde(model=mod1d,statistic=sde.fun1d,R=100,exact=list(Me=1,Mo=1,Va=0.5^2/6),
                 names=c("Me(10)","Mo(10)","Va(10)"))
mc.sde1d
plot(mc.sde1d,index=1)
plot(mc.sde1d,index=2)
plot(mc.sde1d,index=3)

## Example 2 : with Parallel computing
\dontrun{
mod1d <-  snssde1d(drift=f,diffusion=g,x0=5,T=10,M=1000)
## On Windows or Unix
mc.sde1d = MCM.sde(model=mod1d,statistic=sde.fun1d,R=1000,exact=list(Me=1,Mo=1,Va=0.5^2/6),
           names=c("Me(10)","Mo(10)","Va(10)"),parallel="snow",ncpus=parallel::detectCores())
mc.sde1d
## On Unix only
mc.sde1d = MCM.sde(model=mod1d,statistic=sde.fun1d,R=1000,exact=list(Me=1,Mo=1,Va=0.5^2/6),
           names=c("Me(10)","Mo(10)","Va(10)"),parallel="multicore",ncpus=parallel::detectCores())
mc.sde1d
}

## Example 3: (2 dim)
## dX(t) = 1/mu*(theta-X(t)) dt + sqrt(sigma) * dW1(t),
## dY(t) = X(t) dt + 0 * dW2(t)
\dontrun{
## Set the model 2d
mu=0.75;sigma=0.1;theta=2
x0=0;y0=0;init=c(x=0,y=0)
f <- expression(1/mu*(theta-x), x)  
g <- expression(sqrt(sigma),0)
OUI <- snssde2d(drift=f,diffusion=g,M=1000,Dt=0.01,x0=init)

## function of the statistic(s) of interest.
sde.fun2d <- function(data, i){
  d <- data[i,]
  return(c(mean(d$x),mean(d$y),var(d$x),var(d$y),cov(d$x,d$y)))
}
## Monte-Carlo at time = 5
mc.sde2d_a = MCM.sde(model=OUI,statistic=sde.fun2d,R=100,time=5,
                  parallel="snow",ncpus=parallel::detectCores())
mc.sde2d_a
## Monte-Carlo at time = 10
mc.sde2d_b = MCM.sde(model=OUI,statistic=sde.fun2d,R=100,time=10,
                  parallel="snow",ncpus=parallel::detectCores())
mc.sde2d_b

## Compared with exact values at time 5 and 10
E_x <- function(t) theta+(x0-theta)*exp(-t/mu)
V_x <- function(t) 0.5*sigma*mu *(1-exp(-2*(t/mu)))
E_y <- function(t) y0+theta*t+(x0-theta)*mu*(1-exp(-t/mu))
V_y <- function(t) sigma*mu^3*((t/mu)-2*(1-exp(-t/mu))+0.5*(1-exp(-2*(t/mu))))
cov_xy <- function(t) 0.5*sigma*mu^2 *(1-2*exp(-t/mu)+exp(-2*(t/mu)))

## at time=5
mc.sde2d_a = MCM.sde(model=OUI,statistic=sde.fun2d,R=100,time=5,
    exact=list(m1=E_x(5),m2=E_y(5),S1=V_x(5),S2=V_y(5),C12=cov_xy(5)),
	parallel="snow",ncpus=parallel::detectCores())
mc.sde2d_a
plot(mc.sde2d_a,index=1)
plot(mc.sde2d_a,index=2)
## at time=10
mc.sde2d_b = MCM.sde(model=OUI,statistic=sde.fun2d,R=100,time=10,
    exact=list(m1=E_x(10),m2=E_y(10),S1=V_x(10),S2=V_y(10),C12=cov_xy(10)),
	parallel="snow",ncpus=parallel::detectCores())
mc.sde2d_b
plot(mc.sde2d_b,index=1)
plot(mc.sde2d_b,index=2)
}

## Example 4: (3 dim)
## dX(t) = sigma*(Y(t)-X(t)) dt + 0.1 * dW1(t)
## dY(t) = (rho*X(t)-Y(t)-X(t)*Z(t)) dt + 0.1 * dW2(t)
## dZ(t) = (X(t)*Y(t)-bet*Z(t)) dt + 0.1 * dW3(t)
## W1(t), W2(t) and W3(t) are three correlated Brownian motions with Sigma
\dontrun{
## Set the model 3d
sigma=10;rho=28; bet=8/3
f <- expression(sigma*(y-x),rho*x-y-x*z,x*y-bet*z)
g <- expression(0.1,0.1,0.1)
# correlation matrix
Sigma <-matrix(c(1,0.3,0.5,0.3,1,0.2,0.5,0.2,1),nrow=3,ncol=3) 
mod3d <- snssde3d(x0=rep(0,3),drift=f,diffusion=g,M=1000,Dt=0.01,corr=Sigma)

## function of the statistic(s) of interest.
sde.fun3d <- function(data, i){
  d <- data[i,]
  return(c(mean(d$x),mean(d$y),mean(d$z)))
}
## Monte-Carlo at time = 10
mc.sde3d = MCM.sde(mod3d,statistic=sde.fun3d,R=100,parallel="snow",ncpus=parallel::detectCores())
mc.sde3d
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Monte-Carlo}
\keyword{sde}
\keyword{ts}
\keyword{mts}
\keyword{parallel}


