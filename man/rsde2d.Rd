

\name{rsde2d}
\alias{rsde2d}
\alias{rsde2d.default}
\alias{dsde2d}
\alias{dsde2d.default}
\alias{print.dsde2d}
\alias{plot.dsde2d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Approximate transitional densities and random generation for 2-D SDE's}
\description{
Transition density and random generation for the joint and marginal of \code{(X(t-s),Y(t-s) | X(s)=x0,Y(s)=y0)} of the SDE's 2-d.
}
\usage{
rsde2d(object, \dots)
dsde2d(object, ...)

\method{rsde2d}{default}(object, at, \dots)

\method{dsde2d}{default}(object, pdf=c("Joint","Marginal"), at, ...)
\method{plot}{dsde2d}(x,display=c("persp","rgl","image","contour"),hist=FALSE,...)
}
%- maybe also 'usage' for other xs documented here.
\arguments{
  \item{object}{an object inheriting from class \code{\link{snssde2d}} and \code{\link{bridgesde2d}}.}
  \item{at}{time between \code{s=t0} and \code{t=T}. The default \code{at = T}.}
  \item{pdf}{probability density function \code{Joint} or \code{Marginal}.}  
  \item{x}{an object inheriting from class \code{dsde2d}.}
  \item{display}{display plots.} 
  \item{hist}{if \code{hist=TRUE} plot histogram. Based on \code{\link[MASS]{truehist}} function.} 
  \item{\dots}{potentially potentially arguments to be passed to methods, such as \code{\link[stats]{density}} for marginal density and \code{\link[MASS]{kde2d}} fro joint density.}
}
\details{
The function \code{rsde2d} returns a \code{M} random variable \eqn{x_{t=at},y_{t=at}}{x(t=at),y(t=at)} realize at time \eqn{t=at}{t=at}.

\if{html}{\figure{r2d.png}{fig03}}

And \code{dsde2d} returns a bivariate density approximation for \code{(X(t-s),Y(t-s) | X(s)=x0,Y(s)=y0)}. with \eqn{ t=at }{t=at} is a fixed time between \code{t0} and \code{T}.

\if{html}{\figure{tranden2d.png}{fig04}}

An overview of this package, see \code{\link[utils]{browseVignettes}('Sim.DiffProc')} for more informations.
}
\value{
\code{dsde2d} {gives the bivariate density approximation for \code{(X(t-s),Y(t-s) | X(s)=x0,Y(s)=y0)}.}
\code{rsde2d} {generates random of the couple \code{(X(t-s),Y(t-s) | X(s)=x0,Y(s)=y0)}.}
}
\author{
A.C. Guidoum, K. Boukhetala.
}

\seealso{
\code{\link[MASS]{kde2d}} Two-dimensional kernel density estimation in "MASS" package.

\code{\link[ks]{kde}} Kernel density estimate for 1- to 6-dimensional data in "ks" package.

\code{\link[sm]{sm.density}} Nonparametric density estimation in one, two or three dimensions in "sm" package.

\code{\link[yuima]{rng}} random number generators in "yuima" package. 

\code{\link[DiffusionRgqd]{BiGQD.density}} Generate the transition density of a bivariate generalized quadratic diffusion model (2D GQD).
}
\examples{
## Example:1
set.seed(1234)

# SDE's 2d
fx <- expression(3*(2-y),2*x)
gx <- expression(1,y)
mod2d <- snssde2d(drift=fx,diffusion=gx,x0=c(1,2),M=1000)

# random 
r2d <- rsde2d(mod2d,at=0.5)
summary(r2d)

# Marginal density 

denM <- dsde2d(mod2d,pdf="M", at=0.5)
denM
plot(denM)

# Joint density
denJ <- dsde2d(mod2d,pdf="J",n=200, at= 0.5,lims=c(-3,4,0,6))
denJ
plot(denJ)
plot(denJ,display="contour")

## Example 2: Bivariate Transition Density of 2 Brownian motion (W1(t),W2(t)) in [0,1]

\dontrun{
B2d <- snssde2d(drift=rep(expression(0),2),diffusion=rep(expression(1),2),
       M=10000)
for (i in seq(B2d$Dt,B2d$T,by=B2d$Dt)){
plot(dsde2d(B2d, at = i,lims=c(-3,3,-3,3),n=100),
   display="contour",main=paste0('Transition Density \n t = ',i))
}
}

## Example 3: 

\dontrun{
fx <- expression(4*(-1-x)*y , 4*(1-y)*x )
gx <- expression(0.25*y,0.2*x)
mod2d1 <- snssde2d(drift=fx,diffusion=gx,x0=c(x0=1,y0=-1),
      M=5000,type="str")

# Marginal transition density
for (i in seq(mod2d1$Dt,mod2d1$T,by=mod2d1$Dt)){
plot(dsde2d(mod2d1,pdf="M", at = i),main=
      paste0('Marginal Transition Density \n t = ',i))
}

# Bivariate transition density
for (i in seq(mod2d1$Dt,mod2d1$T,by=mod2d1$Dt)){
plot(dsde2d(mod2d1, at = i,lims=c(-1,2,-1,1),n=100),
    display="contour",main=paste0('Transition Density \n t = ',i))
}
}

## Example 4: Bivariate Transition Density of 2 bridge Brownian motion (W1(t),W2(t)) in [0,1]

\dontrun{
B2d <- bridgesde2d(drift=rep(expression(0),2),
  diffusion=rep(expression(1),2),M=5000)
for (i in seq(0.01,0.99,by=B2d$Dt)){ 
plot(dsde2d(B2d, at = i,lims=c(-3,3,-3,3),
 n=100),display="contour",main=
 paste0('Transition Density \n t = ',i))
}
}

## Example 5: Bivariate Transition Density of bridge 
## Ornstein-Uhlenbeck process and its integral in [0,5]
## dX(t) = 4*(-1-X(t)) dt + 0.2 dW1(t)
## dY(t) = X(t) dt + 0 dW2(t)
## x01 = 0 , y01 = 0
## x02 = 0, y02 = 0 
\dontrun{
fx <- expression(4*(-1-x) , x)
gx <- expression(0.2 , 0)
OUI <- bridgesde2d(drift=fx,diffusion=gx,Dt=0.005,M=1000)
for (i in seq(0.01,4.99,by=OUI$Dt)){
plot(dsde2d(OUI, at = i,lims=c(-1.2,0.2,-2.5,0.2),n=100),
 display="contour",main=paste0('Transition Density \n t = ',i))
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{sde}
\keyword{ts}
\keyword{mts}
\keyword{random generators}
\keyword{transition density}
