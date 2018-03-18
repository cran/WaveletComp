\name{FourierRand}
\alias{FourierRand}

\title{Simulation of a surrogate for a given time series x by Fourier randomization}

\description{
It simulates a surrogate with similar Fourier spectrum for a given time series x analyzed 
by wavelet transformation using either function \code{analyze.wavelet} or function \code{analyze.coherency}.
A set of surrogates is used for significance assessment to test the hypothesis of equal periodic components. 

This code has been developed by Huidong Tian (archived R package \code{WaveletCo}).
}
\usage{
FourierRand(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the given time series}
}
% \details{ }
\value{
A surrogate series for x is returned which has the same length and a similar Fourier spectrum. 
}
\references{
Tian, H., and Cazelles, B., 2012. \code{WaveletCo}. 
Available at \url{https://cran.r-project.org/src/contrib/Archive/WaveletCo/}, archived April 2013; accessed July 26, 2013.
}
\author{Code developed by Huidong Tian.}
% \note{ }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{analyze.wavelet}}, \code{\link{analyze.coherency}}, \code{\link{SurrogateData}}, \code{\link{AR}}, \code{\link{ARIMA}}
}
% \examples{ }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
