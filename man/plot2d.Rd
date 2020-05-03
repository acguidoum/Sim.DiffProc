\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\name{plot2d}
\alias{plot2d}
\alias{lines2d}
\alias{points2d}
\alias{plot3D}
\alias{plot2d.default}
\alias{lines2d.default}
\alias{points2d.default}
\alias{plot3D.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plotting for Class SDE}
\description{Generic function for plotting.}
\usage{
\method{plot2d}{default}(x, \dots)
\method{lines2d}{default}(x, \dots)
\method{points2d}{default}(x, \dots)
\method{plot3D}{default}(x, display = c("persp","rgl"), \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object inheriting from class \code{\link{snssde2d}}, \code{\link{snssde3d}}, \code{\link{bridgesde2d}} and \code{\link{bridgesde3d}}.}
  \item{display}{\code{"persp"} perspective or \code{"rgl"} plots.}   
  \item{\dots}{other graphics parameters, see \code{\link[graphics]{par}} in package "graphics", \code{\link[scatterplot3d]{scatterplot3d}} in package "scatterplot3d" and 
  \code{\link[rgl]{plot3d}} in package "rgl".}
}
\details{
The 2 and 3-dim plot of class sde.
}
\author{
A.C. Guidoum, K. Boukhetala.
}
\examples{
## Example 1: 
set.seed(1234)

fx <- rep(expression(0),2)
gx <- rep(expression(1),2)

res <- snssde2d(drift=fx,diffusion=gx,N=5000)
plot2d(res,type="l")

## Example 2: 
set.seed(1234)

fx <- rep(expression(0),3)
gx <- rep(expression(1),3)

res <- snssde3d(drift=fx,diffusion=gx,N=5000)
plot3D(res,display="persp")
plot3D(res,display="rgl")
}
