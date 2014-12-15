\name{wc}
\alias{wc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cross-wavelet transformation, wavelet coherence computation, and a simulation algorithm}
\description{
This function provides Morlet cross-wavelet transformation results of the given two time series, 
performed within the lower-order functions \code{WaveletCoherency} and \code{WaveletTransform} 
subject to criteria concerning the time and frequency resolution, an (optional) lower and/or upper Fourier period, 
and a variety of filtering methods for the coherence computation.
It performs a simulation algorithm to test against a specified alternative hypothesis, which can be chosen from a variety of options,
and provides p-values. The selected model will be fitted to the data and simulated according to estimated parameters in order to provide surrogate time series.
This function is called by function \code{analyze.coherency}.

The name and parts of the layout were inspired by a similar function developed by 
Huidong Tian and Bernard Cazelles (archived R package \code{WaveletCo}). 
The major part of the code for the computation of the cone of influence, and the code for Fourier-randomized surrogate time series have been adopted from Huidong Tian.
The implementation of a choice of filtering windows for the computation of
the wavelet coherence was inspired by Luis Aguiar-Conraria and Maria Joana Soares (\code{GWPackage}).
}
\usage{
wc(x, y, start = 1, dt = 1, dj = 1/20, 
   lowerPeriod = 2*dt, upperPeriod = floor(length(x)/3)*dt, 
   window.type.t = 1, window.type.s = 1, 
   window.size.t = 5, window.size.s = 1/4,
   make.pval = T, method = "white.noise", params = NULL, 
   n.sim = 100, save.sim = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the time series x to be analyzed}
  \item{y}{the time series y to be analyzed (of the same length as x}
  \item{start}{starting value (in time units) at the time axis (for the computation of the cone of influence)}
  \item{dt}{time resolution, i.e. sampling resolution on time domain, 1/dt = number of intervals per time unit. Default: 1.}
  \item{dj}{frequency resolution, i.e. sampling resolution on frequency domain, 1/dj = number of suboctaves (voices per octave). Default: 1/20.}
  \item{lowerPeriod}{lower Fourier period (in time units) for wavelet decomposition. Default: 2*dt.}
  \item{upperPeriod}{upper Fourier period (in time units) for wavelet decomposition. Default: (floor of one third of time series length)*dt.}
  \item{window.type.t}{type of window for smoothing in time direction, select from:  
                 \tabular{rllll}{
                     \tab  0 \tab ("none") \tab : \tab no smoothing in time direction \cr
                     \tab  1 \tab ("bar")  \tab : \tab Bartlett \cr
                     \tab  2 \tab ("tri")  \tab : \tab Triangular (Non-Bartlett) \cr
                     \tab  3 \tab ("box")  \tab : \tab Boxcar (Rectangular, Dirichlet)\cr
                     \tab  4 \tab ("han")  \tab : \tab Hanning \cr
                     \tab  5 \tab ("ham")  \tab : \tab Hamming \cr
                     \tab  6 \tab ("bla")  \tab : \tab Blackman 
                     } 
                  Default: 1 = "bar".
                }
  \item{window.type.s}{type of window for smoothing in scale (period) direction, select from:
                 \tabular{rllll}{
                     \tab  0 \tab ("none") \tab : \tab no smoothing in scale (period) direction \cr
                     \tab  1 \tab ("bar") \tab : \tab Bartlett \cr
                     \tab  2 \tab ("tri") \tab : \tab Triangular (Non-Bartlett) \cr
                     \tab  3 \tab ("box") \tab : \tab Boxcar (Rectangular, Dirichlet)\cr
                     \tab  4 \tab ("han") \tab : \tab Hanning \cr
                     \tab  5 \tab ("ham") \tab : \tab Hamming \cr
                     \tab  6 \tab ("bla") \tab : \tab Blackman 
                     } 
                 Default: 1 = "bar".
                 }
  \item{window.size.t}{size of the window used for smoothing in time direction in units of 1/dt. Default: 5, which together with dt=1 defines a window of length 5*(1/dt) = 5. Windows of even-numbered sizes are extended by 1.}
  \item{window.size.s}{size of the window used for smoothing in scale direction in units of 1/dj. Default: 1/4, which together with dj=1/20 defines a window of length (1/4)*(1/dj) = 5. Windows of even-numbered sizes are extended by 1.}
  \item{make.pval}{Compute p-values? Logical. Default: TRUE.}
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
  
                Default which includes:
  
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
  \item{n.sim}{number of simulations. Default: 100.}
  \item{save.sim}{Shall simulations be saved on the output list? Logical. Default: FALSE.}
}
% \details{ }
\value{A list with the following elements:
  \item{Wave.xy}{(complex-valued) cross-wavelet transform (analogous to Fourier cross-frequency spectrum, and to the covariance in statistics)} 
  \item{Angle}{phase difference, i.e. phase lead of x over y (= phase.x-phase.y)}
  \item{sWave.xy}{smoothed (complex-valued) cross-wavelet transform}
  \item{sAngle}{phase difference, i.e. phase lead of x over y, affected by smoothing} 
  \item{Power.xy}{cross-wavelet power (analogous to Fourier cross-frequency power spectrum)}
  \item{Power.xy.avg}{average cross-wavelet power in the frequency domain (averages over time)} 
  \item{Power.xy.pval}{p-values of cross-wavelet power}
  \item{Power.xy.avg.pval}{p-values of average cross-wavelet power}
  \item{Coherency}{(complex-valued) wavelet coherency of series x over series y in the time/frequency domain, affected by smoothing (analogous to Fourier coherency, and to the coefficient of correlation in statistics)}       
  \item{Coherence}{wavelet coherence (analogous to Fourier coherence, and to the coefficient of determination in statistics (affected by smoothing)}
  \item{Coherence.avg}{average wavelet coherence in the frequency domain (averages across time)}
  \item{Coherence.pval}{p-values of wavelet coherence}
  \item{Coherence.avg.pval}{p-values of average wavelet coherence}
  \item{Wave.x, Wave.y}{(complex-valued) wavelet transforms of series x and y}
  \item{Phase.x, Phase.y}{phases of series x and y}
  \item{Ampl.x, Ampl.y}{amplitudes of series x and y}
  \item{Power.x, Power.y}{wavelet power of series x and y}
  \item{Power.x.avg, Power.y.avg}{average wavelet power of series x and y, averages across time} 
  \item{Power.x.pval, Power.y.pval}{p-values of wavelet power of series x and y}
  \item{Power.x.avg.pval, Power.y.avg.pval}{p-values of average wavelet power of series x and y}
  \item{sPower.x, sPower.y}{smoothed wavelet power of series x and y}
  \item{Period}{the Fourier periods (in time units)}
  \item{Scale}{the scales}
  \item{coi.1, coi.2}{borders of the region where the wavelet transforms are not influenced by edge effects (cone of influence)}
  \item{nc}{number of columns/time steps}
  \item{nr}{number of rows/scales/Fourier periods} 
  \item{axis.1}{tick levels corresponding to time steps}
  \item{axis.2}{tick levels corresponding to Fourier periods (= log2(Period))}
  \item{series.sim}{a data frame of the series simulated as surrogates for the (detrended) time series (if both make.pval = TRUE and save.sim = TRUE.)}
}
\references{

Aguiar-Conraria L., and Soares M.J., 2011. 
Business cycle synchronization and the Euro: A wavelet analysis. 
Journal of Macroeconomics 33 (3), 477--489.

Aguiar-Conraria L., and Soares M.J., 2011. 
The Continuous Wavelet Transform: A Primer. 
NIPE Working Paper Series 16/2011.

Aguiar-Conraria L., and Soares M.J., 2012. \code{GWPackage}. Available at \url{http://sites.google.com/site/aguiarconraria/joanasoares-wavelets}; accessed September 4, 2013.

Cazelles B., Chavez M., Berteaux, D., Menard F., Vik J.O., Jenouvrier S., and Stenseth N.C., 2008.
Wavelet analysis of ecological time series.
Oecologia 156, 287--304.

Liu P.C., 1994. 
Wavelet spectrum analysis and ocean wind waves. 
In: Foufoula-Georgiou E., and Kumar P., (eds.), Wavelets in Geophysics, Academic Press, San Diego, 151--166.

Tian, H., and Cazelles, B., 2012. \code{WaveletCo}. Available at \url{http://cran.r-project.org/src/contrib/Archive/WaveletCo/}, archived April 2013; accessed July 26, 2013.

Torrence C., and Compo G.P., 1998. 
A practical guide to wavelet analysis.
Bulletin of the American Meteorological Society 79 (1), 61--78.

Veleda D., Montagne R., and Araujo M., 2012. 
Cross-Wavelet Bias Corrected by Normalizing Scales. 
Journal of Atmospheric and Oceanic Technology 29, 1401--1408.
}
\author{ Angi Roesch and Harald Schmidbauer; credits are also due to Huidong Tian, Bernard Cazelles, Luis Aguiar-Conraria, and Maria Joana Soares.}
% \note{ }

%% ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{ }
% \examples{ }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }