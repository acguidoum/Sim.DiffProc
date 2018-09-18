\name{NEWS}
\title{\pkg{Sim.DiffProc} News}
\encoding{UTF-8}

\section{CHANGES IN \pkg{Sim.DiffProc} VERSION 4.2}{
    \itemize{
	  \item{Release date: 2018-10-15.}
      \item{Minor bug fixes: R (>= 3.5.0).}
    }
}

\section{CHANGES IN \pkg{Sim.DiffProc} VERSION 4.1}{
    \itemize{
	  \item{Release date: 2018-09-10.}
      \item{Addition a new function \code{TEX.sde()}: converting \pkg{Sim.DiffProc} objects to LaTeX.}
      \item{Many help-file improvements.}
      \item{Manual update.}
    }
}

\section{CHANGES IN \pkg{Sim.DiffProc} VERSION 4.0}{
    \itemize{
	  \item{Release date: 2017-11-28.}
      \item{Manual update.}
      \item{Minor bug fixes: R (>= 3.4.0).}
    }
}

\section{CHANGES IN \pkg{Sim.DiffProc} VERSION 3.8}{
    \itemize{
	  \item{Release date: 2017-09-16.}
      \item{Addition a new function \code{MCM.sde()}: Parallel Monte-Carlo of a statistic applied to SDE's.}
      \item{Addition a new function \code{MEM.sde()}: Moment equations Methods for SDE's.}
      \item{Modified functions: \code{fptsde1d()}, \code{fptsde2d()} and \code{fptsde3d()}, for simulation and approximation the density of first-passage-time.}
      \item{Modified \code{summary()} function for \code{bridgesde1d()}, \code{bridgesde2d()} and \code{bridgesde3d()} classes.}
      \item{Many help-file improvements.}
      \item{Manual update.}
      \item{Minor bug fixes: R (>= 3.4.0).}
    }
}

\section{CHANGES IN \pkg{Sim.DiffProc} VERSION 3.7}{
    \itemize{
	  \item{Release date: 2017-03-25.}
      \item{Addition a new argument in \code{dsde3d()}, for approximation transition density in 3-dim.}
      \item{Minor bug fixes: R (>= 3.3.3).}
    }
}

\section{CHANGES IN \pkg{Sim.DiffProc} VERSION 3.6}{
    \itemize{
	 \item{Release date: 2016-11-13.}
     \item{Minor bug fixes: R (>= 3.3.3).}
    }
}

\section{CHANGES IN \pkg{Sim.DiffProc} VERSION 3.5}{
    \itemize{
	\item{Release date: 2016-11-10.}
      \item{Addition a new functions \code{dsdekd()}, where \code{k=1,2} and \code{3-dim}, for approximation the transition density. }
      \item{Addition a new functions \code{dfptsdekd()}, where \code{k=1,2} and \code{3-dim}, for approximation the density of first-passage-time.	}
      \item{Addition a new functions \code{fptsdekd()}, where \code{k=1,2} and \code{3-dim}, for random generation of first-passage-time.}
      \item{Other minor bug fixes: R (>= 3.3.1).}
    }
}

\section{CHANGES IN \pkg{Sim.DiffProc} VERSION 3.4}{
    \itemize{
	\item{Release date: 2016-10-08.}
      \item{Modified functions \code{rsdekd()}, where \code{k=1,2} and \code{3-dim}.}
      \item{Modified functions \code{fptsdekd()}, where \code{k=1,2} and \code{3-dim}.}
      \item{Other minor bug fixes: R (>= 3.3.1).}
    }
}

\section{CHANGES IN \pkg{Sim.DiffProc} VERSION 3.3}{
    \itemize{
	  \item{Release date: 2016-09-21.}
      \item{Modified \code{summary()} function for \code{snssdekd()}, where \code{k=1,2} and \code{3-dim} classes.}
      \item{Fixed some bugs in \code{"Sim.DiffProc.c"}}
      \item{Other minor bug fixes: R (>= 3.3.0).}
    }
}

\section{CHANGES IN \pkg{Sim.DiffProc} VERSION 3.2}{
    \itemize{
	  \item{Release date: 2016-02-09.}
      \item{Other minor bug fixes: R (>= 3.2.3).}
    }
}

\section{CHANGES IN \pkg{Sim.DiffProc} VERSION 3.1}{
    \itemize{
	  \item{Release date: 2015-12-31.}
      \item{Other minor bug fixes: R (>= 3.2.3).}
    }
}

\section{CHANGES IN \pkg{Sim.DiffProc} VERSION 3.0}{
    \itemize{
	  \item{Release date: 2015-10-30.}
      \item{Other minor bug fixes: R (>= 3.2.1).}
    }
}

\section{CHANGES IN \pkg{Sim.DiffProc} VERSION 2.9}{
    \itemize{
	 \item{Release date: 2014-09-10.}
      \item{Addition a new argument in \code{bridgesdekd()} class, where \code{k=1,2} and \code{3-dim}, possibility of simulating flow paths.}
      \item{Other minor bug fixes: R (>= 3.0.0).}
    }
}

\section{CHANGES IN \pkg{Sim.DiffProc} VERSION 2.8}{
    \itemize{
	  \item{Release date: 2014-08-18.}
      \item{Numerical methods are reprogrammed with \code{C} language 'Quick execution time'.}
      \item{Addition a new argument in \code{snssde2d()} and \code{snssde3d()}, possibility of simulating flow paths.}
      \item{Addition a new argument in \code{fptsdekd()}, where \code{k=1,2} and \code{3-dim}, first-passage-time with time-dependent boundary.  }
      \item{Possibility of displays for \code{fptsde2d()}, \code{fptsde3d()}, \code{rsde2d()} and \code{rsde3d()} functions.}
      \item{Many help-file improvements.}
      \item{Manual update.}
      \item{Other minor bug fixes: R (>= 3.0.0).}
    }
}

\section{CHANGES IN \pkg{Sim.DiffProc} VERSION 2.7}{
    \itemize{
	  \item{Release date: 2014-07-13.}
      \item{Other minor bug fixes: R (>= 3.0.0).}
    }
}

\section{CHANGES IN \pkg{Sim.DiffProc} VERSION 2.6}{
    \itemize{
	  \item{Release date: 2014-05-04.}
      \item{Addition a new type of SDE "Stratonovich" in \code{snssdekd()}, where \code{k=1,2} and \code{3-dim} functions.}
      \item{Addition a new numerice methods "Runge-Kutta Order 1 and 2" add in \code{snssdekd()}, where \code{k=1,2} and \code{3-dim} functions.}
      \item{Addition a new function \code{st.int()}.}
      \item{Addition a new function \code{fitsde()}.}
      \item{Addition a new function \code{fptsdekd()}, where \code{k=1,2} and \code{3-dim}.}
      \item{Addition a new function \code{rsdekd()}, where \code{k=1,2} and \code{3-dim}.}
	  \item{Manual update.}
	  \item{Many help-file improvements.}
	  \item{Other minor bug fixes: R (>= 2.15.1).}
    }
}

\section{CHANGES IN \pkg{Sim.DiffProc} VERSION 2.5}{
    \itemize{
	  \item{Release date: 2012-06-09.}
	  \item{Other minor bug fixes: R (>= 2.15.1).}
    }
}

\section{CHANGES IN \pkg{Sim.DiffProc} VERSION 2.2}{
    \itemize{
	  \item{Release date: 2012-02-13.}
	  \item{Rename the default in \code{"read.xlsx"}, createSheet and removeSheet to \code{sheetName="Sheet1" from "Sheet 1"}.  This makes it consistent with Excel 2007 and 2010 names of an empty workbook.}
    }
}

\section{CHANGES IN \pkg{Sim.DiffProc} VERSION 2.1}{
    \itemize{
	  \item{Release date: 2011-11-18.}
	  \item{A new function \code{Appdcon()} approximated the evolution conditional law of diffusion process by Euler method, Kessler s'methods and Shoji-Ozaki methods.}
    }
}

\section{CHANGES IN \pkg{Sim.DiffProc} VERSION 2.0}{
    \itemize{
	  \item{Release date: 2011-02-13.}
	  \item{Other minor bug fixes: R (>= 2.15.1).}
    }
}

\section{CHANGES IN \pkg{Sim.DiffProc} VERSION 1.0}{
    \itemize{
	  \item{Release date: 2011-01-04.}
	  \item{Version 1.0 released.}
    }
}
