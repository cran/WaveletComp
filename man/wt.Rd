\name{wt}
\alias{wt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Wavelet transformation, and a simulation algorithm}

\description{
This function provides Morlet wavelet transformation results of the given time series, 
performed within the lower-order function \code{WaveletTransform},  
subject to criteria concerning the time and frequency resolution, and an (optional) lower and/or upper Fourier period.
It performs a simulation algorithm to test against a specified alternative hypothesis, which can be chosen from a variety of options,
and provides p-values. The selected model will be fitted to the data and simulated according to estimated parameters in order to provide surrogate time series.
This function is called by function \code{analyze.wavelet}.

The name and parts of the layout were inspired by a similar function developed by Huidong Tian and Bernard Cazelles (archived R package \code{WaveletCo}).
The major part of the code for the computation of the cone of influence, and the code for Fourier-randomized surrogate time series
has been adopted from Huidong Tian.
}
\usage{
wt(x, start = 1, dt = 1, dj = 1/20, 
   lowerPeriod = 2*dt, upperPeriod = floor(length(x)*dt/3), 
   make.pval = T, method = "white.noise", params = NULL, 
   n.sim = 100, save.sim = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the time series to be analyzed}
  \item{start}{starting value (in time units) at the time axis (for the computation of the cone of influence)}
  \item{dt}{time resolution, i.e. sampling resolution on time domain, 1/dt = number of intervals per time step. Default: 1.}
  \item{dj}{frequency resolution, i.e. sampling resolution on frequency domain, 1/dj = number of suboctaves (voices per octave). Default: 1/20.}
  \item{lowerPeriod}{lower Fourier period (in time units) for wavelet decomposition. Default: 2*dt.}
  \item{upperPeriod}{upper Fourier period (in time units) for wavelet decomposition. Default: (floor of one third of time series length)*dt.}
  \item{make.pval}{Compute p-values? Logical. Default: TRUE}
  \item{method}{the method of generating surrogate time series, select from:
                \tabular{rlll}{
                    \tab "white.noise"  \tab : \tab white noise \cr
                    \tab "shuffle"      \tab : \tab shuffling the given time series \cr
                    \tab "Fourier.rand" \tab : \tab time series with a similar spectrum \cr
                    \tab "AR"           \tab : \tab AR(p) \cr  
                    \tab "ARIMA"        \tab : \tab ARIMA(p,0,q) 
                    }
                Default: "white.noise"
                }
  \item{params}{a list of assignments between methods (AR, and ARIMA) and lists of parameter values which apply to surrogates. Default: NULL. 
  
                Default includes:
                
                AR: AR = list(p=1), where:
                
                \tabular{rlll}{
                    \tab p \tab : \tab AR order
                    }
                
                ARIMA: ARIMA = list(p=1, q=1, include.mean=T, sd.fac=1, trim = F, trim.prop = 0.01), where:
                
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
  \item{n.sim}{number of simulations. Default: 100}
  \item{save.sim}{Should simulations be saved on the output list? Logical. Default: FALSE.}
}
% \details{ }
\value{A list with the following elements:
  \item{Wave}{complex wavelet transform of the series}
  \item{Phase}{phases}
  \item{Ampl}{amplitudes}
  \item{Power}{wavelet power in the time/frequency domain}
  \item{Power.avg}{vector of average wavelet power in the frequency domain (averages over time)} 
  \item{Power.pval}{p-values of wavelet power}
  \item{Power.avg.pval}{vector of p-values of average wavelet power}
  \item{Period}{the Fourier periods (in time units)}
  \item{Scale}{the scales}
  \item{coi.1, coi.2}{borders of the region where the wavelet transforms are not influenced by edge effects (cone of influence)}
  \item{nc}{number of columns/time steps}
  \item{nr}{number of rows/scales/Fourier periods}
  \item{axis.1}{tick levels corresponding to time steps}
  \item{axis.2}{tick levels corresponding to Fourier periods (= log2(Period))}
  \item{series.sim}{a data frame of the series which have been simulated as surrogates for the (detrended) time series (if both make.pval = TRUE and save.sim = TRUE.)}
}
\references{

Aguiar-Conraria L., and Soares M.J., 2011. 
The Continuous Wavelet Transform: A Primer. 
NIPE Working Paper Series 16/2011.

Carmona R., Hwang W.-L., and Torresani B., 1998. 
Practical Time Frequency Analysis. Gabor and Wavelet Transforms with an Implementation in S.
Academic Press, San Diego.

Cazelles B., Chavez M., Berteaux, D., Menard F., Vik J.O., Jenouvrier S., and Stenseth N.C., 2008.
Wavelet analysis of ecological time series.
Oecologia 156, 287--304.

Liu Y., Liang X.S., and Weisberg R.H., 2007.
Rectification of the Bias in the Wavelet Power Spectrum. 
Journal of Atmospheric and Oceanic Technology 24, 2093--2102.

Tian, H., and Cazelles, B., 2012. \code{WaveletCo}. Available at \url{http://cran.r-project.org/src/contrib/Archive/WaveletCo/}, archived April 2013; accessed July 26, 2013.

Torrence C., and Compo G.P., 1998. 
A practical guide to wavelet analysis. 
Bulletin of the American Meteorological Society 79 (1), 61--78.
}
\author{Angi Roesch and Harald Schmidbauer; credits are also due to Huidong Tian, and Bernard Cazelles}
% \note{ }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{WaveletTransform}}, \code{\link{analyze.wavelet}} } 
% \examples{ }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }