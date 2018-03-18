\name{AR}
\alias{AR}

\title{Simulation of an AR(p) surrogate for a given time series x}

\description{
It simulates an AR(p) surrogate for the time series x analyzed by wavelet transformation using either function
\code{analyze.wavelet} or function \code{analyze.coherency}. A set of surrogates is used for significance assessment 
to test the hypothesis of equal periodic components.

Simulation is subject to coefficient estimates referring to the given AR order.
}
\usage{
AR(x, params = list(AR = list(p = 1)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the given time series}
  %%%%%%%%%%%%%%%%%
  \item{params}{a list containing an assignment of (the only) parameter value to AR. Default: \code{NULL}. 
  
  Default includes:
  
 \code{AR = list(...)}, a list containing one single element:
        \tabular{rlll}{
        \tab \code{p} \tab : \tab AR order.\cr
        \tab          \tab   \tab Default: \code{1}.
        }
  }              
}
\details{
Basics of model estimation and simulation are provided by application of the R functions \code{pcf} and \code{arima.sim}.
}
\value{
A surrogate series for x is returned which has the same length and properties according to estimates resulting from the given AR order.
}
% \references{ }
\author{Angi Roesch and Harald Schmidbauer}
% \note{ }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{analyze.wavelet}}, \code{\link{analyze.coherency}}, \code{\link{SurrogateData}}, \code{\link{ARIMA}}, \code{\link{FourierRand}}
}
% \examples{ }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }