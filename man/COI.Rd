\name{COI}
\alias{COI}

\title{Computation of the cone of influence}

\description{
Wavelet transformation of a time series includes padding with zeros at the edges.
This function determines the spectrum borders influenced by edge effects.

This code has been developed by Huidong Tian (archived R package \code{WaveletCo}).
}
\usage{
COI(start = start, dt = dt, nc = nc, nr = nr, Period = Period)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{start}{starting value (in time units) at the time axis}
  \item{dt}{time resolution, i.e. sampling resolution on time domain}
  \item{nc}{number of columns/time steps}
  \item{nr}{number of rows/scales}
  \item{Period}{the Fourier periods (in time units)}
}
% \details{ }
\value{a list with the following elements:
       \tabular{rlll}{
       \tab x, y \tab : \tab coordinates of the cone of influence \cr
       \tab axis.1, axis.2 \tab : \tab the axes for plotting purpose
       }
}
\references{
Tian, H., and Cazelles, B., 2012. \code{WaveletCo}. Available at \url{http://cran.r-project.org/src/contrib/Archive/WaveletCo/}, archived April 2013; accessed July 26, 2013.
}
\author{Code developed by Huidong Tian.}
% \note{ }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{analyze.wavelet}}, \code{\link{analyze.coherency}}, \code{\link{wt.image}}, \code{\link{wc.image}}
}
% \examples{ }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }