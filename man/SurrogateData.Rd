\name{SurrogateData}
\alias{SurrogateData}

\title{Simulation of surrogates for a given time series x, subject to the specified method and parameters}

\description{
It simulates a surrogate for the time series x to be analyzed by wavelet transformation using either function 
\code{analyze.wavelet} or function \code{analyze.coherency}. A set of surrogates is used for significance assessment 
to test the hypothesis of equal periodic components. 

Simulation is subject to model/method specification and parameter setting:
Currently, one can choose from a variety of 6 methods (white noise, series shuffling, Fourier randomization, AR, and ARIMA) 
with respective lists of parameters to set.

The name and layout were inspired by a similar function developed by Huidong Tian (archived R package \code{WaveletCo}).
}
\usage{
SurrogateData(x, method = "white.noise", params = list( 
         AR    = list(p = 1), 
         ARIMA = list(p = 1, q = 1, include.mean = TRUE, sd.fac = 1, 
                      trim = FALSE, trim.prop = 0.01)))
}

\arguments{
  \item{x}{the given time series}
  %%%%%%%%%%%%%%%%% 
  \item{method}{the method of generating surrogate time series; select from:
                \tabular{rlll}{
                    \tab \code{"white.noise"}  \tab : \tab white noise \cr
                    \tab \code{"shuffle"}      \tab : \tab shuffling the given time series \cr
                    \tab \code{"Fourier.rand"} \tab : \tab time series with a similar spectrum \cr
                    \tab \code{"AR"}           \tab : \tab AR(p) \cr  
                    \tab \code{"ARIMA"}        \tab : \tab ARIMA(p,0,q) 
                    }
  Default: \code{"white.noise"}.
  }
  \item{params}{a list of assignments between methods (AR, and ARIMA) and lists of parameter values 
  applying to surrogates. Default: \code{NULL}. 
  
  Default includes:
    
  \code{AR = list(p = 1)},\cr where:
  \tabular{rlll}{
  \tab \code{p} \tab : \tab AR order
                }
  \code{ARIMA = list(p = 1, q = 1, include.mean = TRUE, sd.fac = 1,}\cr 
  \code{trim = FALSE, trim.prop = 0.01)},\cr where:                
  \tabular{rlll}{
  \tab \code{p}            \tab : \tab  AR order \cr
  \tab \code{q}            \tab : \tab  MA order \cr
  \tab \code{include.mean} \tab : \tab  Include a mean/intercept term? \cr
  \tab \code{sd.fac}       \tab : \tab  magnification factor to boost the \cr 
  \tab                     \tab   \tab  residual standard deviation \cr  
  \tab \code{trim}         \tab : \tab  Simulate trimmed data? \cr
  \tab \code{trim.prop}    \tab : \tab  high/low trimming proportion
                }
  }      
}
% \details{ }
\value{
A surrogate series for x is returned which has the same length and properties according to estimates 
resulting from the model/method specification and parameter setting. 
}
\references{
Tian, H., and Cazelles, B., 2012. \code{WaveletCo}. 
Available at \url{https://cran.r-project.org/src/contrib/Archive/WaveletCo/}, archived April 2013; accessed July 26, 2013.
}
\author{Angi Roesch and Harald Schmidbauer; credits are also due to Huidong Tian.}
% \note{ }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{analyze.wavelet}}, \code{\link{analyze.coherency}}, \code{\link{AR}}, \code{\link{ARIMA}}, \code{\link{FourierRand}}
}
% \examples{ }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }