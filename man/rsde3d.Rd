

\name{rsde3d}
\alias{rsde3d}
\alias{rsde3d.default}
\alias{dsde3d}
\alias{dsde3d.default}
\alias{print.dsde3d}
\alias{plot.dsde3d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Approximate transitional densities and random generation for 3-D SDE's}
\description{
Transition density and random generation for the joint and marginal of \code{(X(t-s),Y(t-s),Z(t-s) | X(s)=x0,Y(s)=y0,Z(s)=z0)} of the SDE's 3-d.
}
\usage{
rsde3d(object, \dots)
dsde3d(object, ...)

\method{rsde3d}{default}(object, at, \dots)

\method{dsde3d}{default}(object, pdf=c("Joint","Marginal"), at, ...)
\method{plot}{dsde3d}(x,display="rgl",hist=FALSE,...)
}
%- maybe also 'usage' for other xs documented here.
\arguments{
  \item{object}{an object inheriting from class \code{\link{snssde3d}} and \code{\link{bridgesde3d}}.}
  \item{at}{time between \code{s=t0} and \code{t=T}. The default \code{at = T}.}
  \item{pdf}{probability density function \code{Joint} or \code{Marginal}.}  
  \item{x}{an object inheriting from class \code{dsde3d}.}
  \item{display}{display plots.} 
  \item{hist}{if \code{hist=TRUE} plot histogram. Based on \code{\link[MASS]{truehist}} function.} 
  \item{\dots}{potentially arguments to be passed to methods, such as \code{\link[stats]{density}} for marginal density  and \code{\link[sm]{sm.density}} for joint density.}
}
\details{
The function \code{rsde3d} returns a \code{M} random variable \eqn{x_{t=at},y_{t=at},z_{t=at}}{x(t=at),y(t=at),z(t=at)} realize at time \eqn{t=at}{t=at}.

\if{html}{\figure{r3d.png}{fig05}}

And \code{dsde3d} returns a trivariate kernel density approximation for \code{(X(t-s),Y(t-s),Z(t-s) | X(s)=x0,Y(s)=y0,Z(s)=z0)}. with \eqn{t= at }{t=at} is a fixed time between \code{t0} and \code{T}.

\if{html}{\figure{tranden3d.png}{fig06}}

An overview of this package, see \code{\link[utils]{browseVignettes}('Sim.DiffProc')} for more informations.
}
\value{
\code{dsde3d} {gives the trivariate density approximation \code{(X(t-s),Y(t-s),Z(t-s) | X(s)=x0,Y(s)=y0,Z(s)=z0)}.}
\code{rsde3d} {generates random of the \code{(X(t-s),Y(t-s),Z(t-s) | X(s)=x0,Y(s)=y0,Z(s)=z0)}.}
}
\author{
A.C. Guidoum, K. Boukhetala.
}

\seealso{
\code{kde} Kernel density estimate for 1- to 6-dimensional data in "ks" package.

\code{sm.density} Nonparametric density estimation in one, two or three dimensions in "sm" package.

\code{kde3d} Compute a three dimension kernel density estimate in "misc3d" package.

\code{rng} random number generators in "yuima" package. 

\code{rcBS}, \code{rcCIR}, \code{rcOU} and \code{rsOU} in package "sde". 
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
mod3d1 <- snssde3d(x0=c(0,2,-2),drift=fx,diffusion=gx,M=1000,Dt=0.003)

# random at t= 0.75
r3d1 <- rsde3d(mod3d1,at=0.75)
summary(r3d1)

# Marginal transition density at t=0.75, t0=0

denM <- dsde3d(mod3d1,pdf="M",at=0.75)
denM
plot(denM)

# for Joint transition density at t=0.75;t0=0 
# Multiple isosurfaces
\dontrun{
denJ <- dsde3d(mod3d1,pdf="J", at= 0.75)
denJ
plot(denJ,display="rgl")
}

## Example 2: Stratonovich sde
## dX(t) = Y(t)* dt + X(t) o dW1(t)          
## dY(t) = (4*( 1-X(t)^2 )* Y(t) - X(t))* dt + 0.2 o dW2(t)
## dZ(t) = (4*( 1-X(t)^2 )* Z(t) - X(t))* dt + 0.2 o dW3(t)
set.seed(1234)

fx <- expression( y , (4*( 1-x^2 )* y - x), (4*( 1-x^2 )* z - x))
gx <- expression( x , 0.2, 0.2)
mod3d2 <- snssde3d(drift=fx,diffusion=gx,M=1000,type="str")

# random 
r3d2 <- rsde3d(mod3d2)
summary(r3d2)

# Marginal transition density at t=1, t0=0

denM <- dsde3d(mod3d2,pdf="M")
denM
plot(denM)

# for Joint transition density at t=1;t0=0
# Multiple isosurfaces
\dontrun{
denJ <- dsde3d(mod3d2,pdf="J")
denJ
plot(denJ,display="rgl")
}

## Example 3: Tivariate Transition Density of 3 Brownian motion (W1(t),W2(t),W3(t)) in [0,1]

\dontrun{
B3d <- snssde3d(drift=rep(expression(0),3),diffusion=rep(expression(1),3),M=500)
for (i in seq(B3d$Dt,B3d$T,by=B3d$Dt)){
plot(dsde3d(B3d, at = i,pdf="J"),box=F,main=paste0('Transition Density t = ',i))
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
