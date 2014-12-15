\name{ARIMA}
\alias{ARIMA}

\title{Simulation of an ARIMA(p,0,q) surrogate for a given time series x}

\description{
It simulates an ARIMA(p,0,q) surrogate for the time series x analyzed by wavelet transformation using either function
\code{analyze.wavelet} or function \code{analyze.coherency}. A set of surrogates is used for significance assessment to test the hypothesis of equal periodic components. 

Simulation is subject to parameters concerning:
AR and MA order, an indicator of mean/intercept term inclusion, a boosting factor which is applied to the residual standard deviation, an indicator to trim the data according to the given trimming proportion (i.e. to restrict them to the interval between the lower and upper quantile corresponding to the given probability).
}
\usage{
ARIMA(x, params = 
      list(ARIMA = list(p = 1, q = 1, 
                        include.mean = T, sd.fac = 1, 
                        trim = F, trim.prop = 0.01)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the given time series}
  \item{params}{a list containing an assignment of parameter values to ARIMA. Default: NULL. 
  
                Default includes:
                
                ARIMA = list(p=1, q=1, include.mean=T, sd.fac=1, trim = F, trim.prop = 0.01), where:
                
                \tabular{rlll}{  
                    \tab  p            \tab : \tab  AR order \cr
                    \tab  q            \tab : \tab  MA order \cr
                    \tab  include.mean \tab : \tab  Include a mean/intercept term? \cr
                    \tab  sd.fac       \tab : \tab  magnification factor to boost the \cr 
                    \tab               \tab   \tab  residual standard deviation \cr  
                    \tab  trim         \tab : \tab  Simulate trimmed data? \cr
                    \tab  trim.prop    \tab : \tab  high/low trimming proportion
                    }

               }
}
\details{
Basics of model estimation and simulation are provided by application of the R functions \code{arima} and \code{arima.sim}.
}
\value{
A surrogate series for x is returned which has the same length and properties according to estimates resulting from the given ARIMA setting. 
}
% \references{ }
\author{Angi Roesch and Harald Schmidbauer}
% \note{ }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{analyze.wavelet}}, \code{\link{analyze.coherency}}, \code{\link{SurrogateData}}, \code{\link{AR}}, \code{\link{FourierRand}}
}
% \examples{ }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }