\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\name{moment}
\alias{bconfint}
\alias{kurtosis}
\alias{moment}
\alias{Median}
\alias{Mode}
\alias{cv}
\alias{skewness}
\alias{bconfint.default}
\alias{kurtosis.default}
\alias{moment.default}
\alias{cv.default}
\alias{Median.default}
\alias{Mode.default}
\alias{skewness.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Monte-Carlo statistics of SDE's}
\description{Generic function for compute the kurtosis, skewness, median, mode and coefficient of variation (relative variability), moment and confidence interval of class \code{"sde"}.}
\usage{
\method{bconfint}{default}(x, level = 0.95, \dots)
\method{kurtosis}{default}(x, \dots)
\method{moment}{default}(x, order = 1,center = TRUE, \dots)
\method{cv}{default}(x, \dots)
\method{skewness}{default}(x, \dots)
\method{Median}{default}(x, \dots)
\method{Mode}{default}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object inheriting from class \code{"sde"}.}  
  \item{order}{order of moment.}
  \item{center}{if \code{TRUE} is a central moment.}
  \item{level}{the confidence level required.} 
  \item{\dots}{potentially further arguments for (non-default) methods.}
}
\author{
A.C. Guidoum, K. Boukhetala.
}

\examples{
## Example 1:
## dX(t) = 2*(3-X(t)) *dt + dW(t)
set.seed(1234)

f <- expression( 2*(3-x) )
g <- expression( 1 )
mod <- snssde1d(drift=f,diffusion=g,M=10000,T=5)
## Monte-Carlo statistics of 5000 trajectory of X(t) at final time T of 'mod'
summary(mod)
kurtosis(mod)
skewness(mod)
mean(mod)
Median(mod)
Mode(mod)
moment(mod,order=4)
cv(mod)
bconfint(mod,level = 0.95) ## of mean
}
