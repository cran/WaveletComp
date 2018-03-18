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
  \item{start}{starting point in time. Default: \code{start = 1}.}
  %%%%%%%%%%%%%%%%%
  \item{dt}{time resolution, i.e. sampling resolution in the time domain, \code{1/dt} = number of observations per time unit}
  %%%%%%%%%%%%%%%%%
  \item{nc}{number of columns = number of observations = number of observation epochs; "epoch" meaning point in time}
  \item{nr}{number of rows = number of scales (Fourier periods)}
  %%%%%%%%%%%%%%%%%
  \item{Period}{the Fourier periods (measured in time units determined by \code{dt})}
}
% \details{ }
\value{a list with the following elements:
       \tabular{rlll}{
       \tab \code{x}, \code{y}   \tab : \tab coordinates of the cone of influence \cr
       \tab \code{axis.1}        \tab : \tab tick levels corresponding to the time steps \cr
       \tab                      \tab   \tab used for wavelet transformation: \cr
       \tab                      \tab   \tab \code{1, 1+dt, 1+2dt, ...} if \code{start = 1} \cr
       \tab \code{axis.2}        \tab : \tab tick levels corresponding to the log of Fourier periods: \cr 
       \tab                      \tab   \tab \code{log2(Period)}
       }
}
\references{
Tian, H., and Cazelles, B., 2012. \code{WaveletCo}. 
Available at \url{https://cran.r-project.org/src/contrib/Archive/WaveletCo/}, archived April 2013; accessed July 26, 2013.
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