\name{WaveletTransform}
\alias{WaveletTransform}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Morlet wavelet transformation of a time series}
\description{ 
It computes the Morlet wavelet transformation of a given time series, subject to criteria concerning:
the time and frequency resolution, an (optional) lower and/or upper Fourier period.

The output is further processed by higher-order functions \code{wt}, \code{WaveletCoherency} and 
\code{wc}, and can be retrieved from \code{analyze.wavelet} and \code{analyze.coherency}.

The name and layout were inspired by a similar function developed by Huidong Tian and Bernard Cazelles (archived R package \code{WaveletCo}).
}
\usage{
WaveletTransform(x, dt = 1, dj = 1/20, 
                 lowerPeriod = 2*dt, upperPeriod = floor(length(x)*dt/3))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the time series to be analyzed}
  \item{dt}{time resolution, i.e. sampling resolution on time domain. Default: 1.}
  \item{dj}{frequency resolution, i.e. sampling resolution on frequency domain, 1/dj = number of voices per octave. Default: 1/20.}
  \item{lowerPeriod}{lower Fourier period (in time units) for wavelet decomposition. Default: 2*dt.}
  \item{upperPeriod}{upper Fourier period (in time units) for wavelet decomposition. Default: (floor of one third of time series length)*dt.}
}
% \details{ }
\value{A list of class \code{analyze.wavelet} with the following elements:
  \item{Wave}{complex wavelet transform of the series}
  \item{Phase}{phases}
  \item{Ampl}{amplitudes}
  \item{Period}{the Fourier periods (in time units)}
  \item{Scale}{the scales}
  \item{Power}{wavelet power in the time/frequency domain}
  \item{nc}{number of columns/time steps}
  \item{nr}{number of rows/scales/periods}
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
\author{Angi Roesch and Harald Schmidbauer; credits are also due to Huidong Tian, and Bernard Cazelles.}
% \note{ }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{wt}}, \code{\link{WaveletCoherency}}, \code{\link{wc}}, \code{\link{analyze.wavelet}}, \code{\link{analyze.coherency}}
}
% \examples{ }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }