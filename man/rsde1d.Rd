\name{rsde1d}
\alias{rsde1d}
\alias{rsde1d.default}
\alias{dsde1d}
\alias{dsde1d.default}
\alias{print.dsde1d}
\alias{plot.dsde1d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Approximate transitional densities and random generation for 1-D SDE
}
\description{
Transition density and random generation for \code{X(t-s) | X(s)=x0} of the 1-dim SDE.
}
\usage{
rsde1d(object, ...)
dsde1d(object, ...)

\method{rsde1d}{default}(object, at, ...)

\method{dsde1d}{default}(object, at, ...)
\method{plot}{dsde1d}(x,hist=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object inheriting from class \code{\link{snssde1d}} and \code{\link{bridgesde1d}}.
}
  \item{at}{
time between \code{s=t0} and \code{t=T}. The default \code{at = T}.
}
  \item{x}{
an object inheriting from class \code{dsde1d}.
}
  \item{hist}{if \code{hist=TRUE} plot histogram. Based on \code{\link[MASS]{truehist}} function.} 
  \item{\dots}{
potentially arguments to be passed to methods, such as \code{\link[stats]{density}} for kernel density.
}
}
\details{
The function \code{rsde1d} returns a \code{M} random variable \eqn{x_{t=at}}{x(t=at)} realize at time \eqn{t=at}{t=at} defined by : 
\deqn{x_{ t=at } = \{ t \geq 0 ; x = X_{ t=at } \}}{x(t=at) = { t>=0 ; x = X(t=at) }}

\if{html}{\figure{r1d.png}{fig01}}

And \code{dsde1d} returns a transition density approximation for \code{X(t-s) | X(s)=x0}.
with \eqn{t= at }{t=at} is a fixed time between \code{t0} and \code{T}.

\if{html}{\figure{tranden1d.png}{fig02}}

An overview of this package, see \code{\link[utils]{browseVignettes}('Sim.DiffProc')} for more informations.
}
\value{
\code{dsde1d} {gives the transition density estimate of \code{X(t-s) | X(s)=x0}.}
\code{rsde1d} {generates random of \code{X(t-s) | X(s)=x0}.}
}
\author{
A.C. Guidoum, K. Boukhetala.
}
\seealso{
\code{\link[stats]{density}} Kernel density estimation in "stats" package. 

\code{kde} Kernel density estimate for 1- to 6-dimensional data in "ks" package.

\code{sm.density} Nonparametric density estimation in one, two or three dimensions in "sm" package.

\code{rng} random number generators in "yuima" package. 

\code{dcSim} Pedersen's simulated transition density in "sde" package.

\code{rcBS}, \code{rcCIR}, \code{rcOU} and \code{rsOU}  in package "sde".

\code{dcBS}, \code{dcCIR}, \code{dcOU} and \code{dsOU}  in package "sde".

\code{GQD.density} Generate the transition density of a scalar generalized quadratic diffusion in "DiffusionRgqd" package.
 
}
\examples{

## Example 1:  
## dX(t) = (-2*(X(t)<=0)+2*(X(t)>=0)) *dt + 0.5 * dW(t)
set.seed(1234)

f <- expression(-2*(x<=0)+2*(x>=0))
g <- expression(0.5)
res1 <- snssde1d(drift=f,diffusion=g,M=5000)
x <- rsde1d(res1, at = 1)
summary(x)
dens1 <-  dsde1d(res1, at = 1)
dens1
plot(dens1,main="Transition density of X(t=1)|X(s=0)=0") # kernel estimated
plot(dens1,hist=TRUE) # histogramme

## Example 2:
## Transition density of standard Brownian motion W(t) at time = 0.5
set.seed(1234)

f <- expression(0)
g <- expression(1)
res2 <- snssde1d(drift=f,diffusion=g,M=5000)
plot(dsde1d(res2, at = 0.5),dens=function(x) dnorm(x,0,sqrt(0.5)))
plot(dsde1d(res2, at = 0.5),dens=function(x) dnorm(x,0,sqrt(0.5)),hist=TRUE)

## Example 3: Transition density of Brownian motion W(t) in [0,1]

\dontrun{
for (i in seq(res2$t0,res2$T,by=res2$Dt)){
plot(dsde1d(res2, at = i),main=paste0('Transition Density \n t = ',i))
}
}

## Example 4:
## Transition density of bridge Brownian motion W(t) at time = 0.25 and 0.75 
set.seed(1234)
\dontrun{
f <- expression(0)
g <- expression(1)
Bd <- bridgesde1d(drift=f,diffusion=g,M=5000)
Bd
plot(dsde1d(Bd, at = 0.25))         ## Transition Density at time=0.25
plot(dsde1d(Bd, at = 0.75),add=TRUE)## Transition Density at time=0.75
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{sde}
\keyword{ts}
\keyword{mts}
\keyword{random generators}
\keyword{transition density}
