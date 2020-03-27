

\name{Irates}
\alias{Irates}
\docType{data}
\title{
Monthly Interest Rates
}
\description{
  monthly observations from 1946--12 to 1991--02

 \emph{number of observations} :  531

 \emph{observation} :  country

 \emph{country} :  United--States
}
\usage{data(Irates)}
\format{A time serie containing :
  \describe{
    \item{r1}{interest rate for a maturity of 1 months  (\% per year). }
    \item{r2}{interest rate for a maturity of 2 months  (\% per year). }
    \item{r3}{interest rate for a maturity of 3 months  (\% per year). }
    \item{r5}{interest rate for a maturity of 5 months  (\% per year). }
    \item{r6}{interest rate for a maturity of 6 months  (\% per year). }
    \item{r11}{interest rate for a maturity of 11 months  (\% per year). }
    \item{r12}{interest rate for a maturity of 12 months  (\% per year). }
    \item{r36}{interest rate for a maturity of 36 months  (\% per year). }
    \item{r60}{interest rate for a maturity of 60 months  (\% per year). }
    \item{r120}{interest rate for a maturity of 120 months  (\% per year). }
  }
}
\source{
McCulloch, J.H. and Kwon, H.C. (1993). U.S. term structure data, 1947--1991, Ohio State Working Paper 93--6, Ohio State University, Columbus

These datasets \code{Irates} are in package "Ecdat".
}

\references{
Croissant, Y. (2014). Ecdat: Data sets for econometrics. R package version 0.2--5.
}

\examples{
data(Irates)
rates <- Irates[,"r1"]
rates <- window(rates, start=1964.471, end=1989.333)

## CKLS modele vs CIR modele 
## CKLS :  dX(t) = (theta1+theta2* X(t))* dt + theta3 * X(t)^theta4 * dW(t)

fx <- expression(theta[1]+theta[2]*x)
gx <- expression(theta[3]*x^theta[4])
fitmod <- fitsde(rates,drift=fx,diffusion=gx,pmle="euler",start = list(theta1=1,theta2=1,
                  theta3=1,theta4=1),optim.method = "L-BFGS-B")
theta <- coef(fitmod)
	
N <- length(rates)
res <- snssde1d(drift=fx,diffusion=gx,M=1000,t0=time(rates)[1],T=time(rates)[N],
                Dt=deltat(rates),x0=rates[1],N=N)

plot(res,type="n",ylim=c(0,35))
lines(rates,col=2,lwd=2)
lines(time(res),apply(res$X,1,mean),col=3,lwd=2)
lines(time(res),apply(res$X,1,bconfint,level=0.95)[1,],col=4,lwd=2)
lines(time(res),apply(res$X,1,bconfint,level=0.95)[2,],col=4,lwd=2)
legend("topleft",c("real data","mean path",
  paste("bound of", 95," confidence")),inset = .01,
  col=2:4,lwd=2,cex=0.8)
}
\keyword{datasets}
