\name{WaveletCoherency}
\alias{WaveletCoherency}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computation of the cross-wavelet power and wavelet coherence of two time series}
\description{
Given two time series x and y having the same length and timestamp, 
this function computes the cross-wavelet power and wavelet coherence applying the Morlet
wavelet, subject to criteria concerning:
the time and frequency resolution, an (optional) lower and/or upper Fourier period, and 
filtering method for the coherence computation.

The output is further processed by the higher-order function \code{wc}
and can be retrieved from \code{analyze.coherency}.

The name and layout were inspired by a similar function developed by Huidong Tian and Bernard Cazelles (archived R package \code{WaveletCo}).
The implementation of a choice of filtering windows for the computation of
the wavelet coherence was inspired by Luis Aguiar-Conraria and Maria Joana Soares (\code{GWPackage}).
}
\usage{
WaveletCoherency(x, y, dt = 1, dj = 1/20, 
                 lowerPeriod = 2*dt, upperPeriod = floor(length(x)*dt/3), 
                 window.type.t = 1, window.type.s = 1, 
                 window.size.t = 5, window.size.s = 1/4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the time series x to be analyzed}
  \item{y}{the time series y to be analyzed (of the same length as x}
  \item{dt}{time resolution, i.e. sampling resolution on time domain, 1/dt = number of intervals per time step. Default: 1.}
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
}
% \details{ }
\value{A list with the following elements:
  \item{Wave.xy}{(complex-valued) cross-wavelet transform (analogous to Fourier cross-frequency spectrum, and to the covariance in statistics)} 
  \item{sWave.xy}{smoothed (complex-valued) cross-wavelet transform}
  \item{Power.xy}{cross-wavelet power (analogous to Fourier cross-frequency power spectrum)} 
  \item{Coherency}{(complex-valued) wavelet coherency of series x over series y in the time/frequency domain, affected by smoothing (analogous to Fourier coherency, and to the coefficient of correlation in statistics)}       
  \item{Coherence}{wavelet coherence (analogous to Fourier coherence, and to the coefficient of determination in statistics (affected by smoothing)}
  \item{Wave.x, Wave.y}{(complex-valued) wavelet transforms of series x and y}
  \item{Phase.x, Phase.y}{phases of series x and y}
  \item{Ampl.x, Ampl.y}{amplitudes of series x and y}
  \item{Power.x, Power.y}{wavelet power of series x and y}
  \item{sPower.x, sPower.y}{smoothed wavelet power of series x and y}
  \item{Period}{the Fourier periods (in time units)}
  \item{Scale}{the scales}
  \item{nc}{number of columns/time steps}
  \item{nr}{number of rows/scales}
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
\author{Angi Roesch and Harald Schmidbauer; credits are also due to Huidong Tian, Bernard Cazelles, Luis Aguiar-Conraria, and Maria Joana Soares.}
% \note{ }
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
   \code{\link{WaveletTransform}}, \code{\link{wc}}, \code{\link{analyze.coherency}}
}
% \examples{ }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }