\name{periodic.series}
\alias{periodic.series}

\title{
Computation of a (deterministic) periodic time series of linearly changing
period.
}

\description{
It computes and returns a sinusoid of a specified length, which has the given initial phase, and linearly changing periods (if requested) starting from a given period length through the given length at the end. There is an option to plot the time series. 
}
\usage{
  periodic.series(start.period = 100, end.period = start.period, 
                  phase = 0, 
                  length = 600, 
                  make.plot = F)
}

\arguments{
  \item{start.period}{period length at start (in steps of time). Default: 100}
  \item{end.period}{period length at end (in steps of time). Default: 50}
  \item{phase}{phase difference (in steps of time), i.e. part of period length which has elapsed relative to the origin. Default: 0}
  \item{length}{number of time steps. Default: 600}
  \item{make.plot}{Plot time series? Logical. Default: FALSE}
}
\details{
This function can be used for illustrating methods and functions. Producing a sinusoid, periodic.series will work best if start.period (and end.period, if different from start.period) are not too small.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
the series as vector
}
% \references{ }
\author{Angi Roesch and Harald Schmidbauer}
% \note{ }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{analyze.wavelet}}, \code{\link{wt.image}}, \code{\link{reconstruct}}
}
\examples{
# The following time series involves periods from 100 through 50:
x = periodic.series(start.period = 100, end.period = 50, make.plot = TRUE)
title("time series with linearly changing period from 100 through 50")

# The following three time series involve three different types of period evolution 
# starting from period 100:
x1 = 0.8*periodic.series(start.period = 100, end.period = 95, phase = 0, length = 1000)
x2 =     periodic.series(start.period = 100, end.period = 100, phase = 0, length = 1000)
x3 = 1.2*periodic.series(start.period = 100, end.period = 105, phase = 0, length = 1000)

ts.plot(x2, ylim = c(-2, +2), xlab = "time", ylab = "series with variable period")
lines(x1, col = "blue")
lines(x3, col = "red")
legend("topleft", legend = c("speeding up (end period = 95)", "period = 100", 
   "slowing down (end period = 105)"), 
   lty = 1, col = c("blue", "black", "red"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }