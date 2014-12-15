\name{WaveletComp-package}
\alias{WaveletComp-package}
\alias{WaveletComp}
\docType{package}
\title{Computational Wavelet Analysis}

\description{Wavelet analysis and reconstruction of time series, cross-wavelets and phase-difference (with filtering options), significance with bootstrap algorithms.}

\details{
\tabular{ll}{
Package: \tab WaveletComp\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-12-15\cr
License: \tab GPL-2\cr
URL: \tab \url{http://www.hs-stat.com/projects/WaveletComp/WaveletComp_guided_tour.pdf} \cr
}

Periodic phenomena of a single time series can be analyzed with function \code{analyze.wavelet}. Results of the analysis (a time/period image of the wavelet power spectrum, plots of the average power, and phase plots for selected periods and a time/period image of phases) can be accessed through various plot functions (\code{wt.image}, \code{wt.avg}, \code{wt.sel.phases}, \code{wt.phase.image}). Function \code{reconstruct} returns the reconstructed time series where reconstruction is according to constraints on significance, period specification, and cone of influence.

The cross-wavelet spectrum and coherency spectrum of two time series can be analyzed with function \code{analyze.coherency}. Results (a time/period image of cross-wavelet power or coherency, plots of averages, plots of phases and phase-differences for selected periods and the time/period image of phase-differences) can be accessed through corresponding functions (\code{wc.image}, \code{wc.avg}, \code{wc.sel.phases}, \code{wc.phasediff.image}).

Detrending of the time series at hand is offered as an option. Wavelet transformations are computed using the Morlet wavelet. Smoothing filters are provided in the case of cross-wavelet transformation to compute wavelet coherency.

Significance is assessed with simulation algorithms, a variety of alternative hypotheses to test is available, for which surrogate time series are provided: white noise, shuffling the given time series, time series with a similar spectrum, AR, and ARIMA.

Names and parts of the layout of some routines were inspired by similar functions developed by 
Huidong Tian and Bernard Cazelles (archived R package \code{WaveletCo}). The basic concept of the simulation algorithm, and of ridge determination build on ideas developed by these authors. The major part of the code for the computation of the cone of influence, and the code for Fourier-randomized surrogate time series
has been adopted from Huidong Tian. 
The implementation of a choice of filtering windows for the computation of the wavelet coherence was inspired by Luis Aguiar-Conraria and Maria Joana Soares (\code{GWPackage}).

Cross-wavelet and coherence computation, the simulation algorithm, and ridge determination build heavily on the use of matrices in order to minimize computation time in R. 
}
\author{Angi Roesch and Harald Schmidbauer; credits are also due to Huidong Tian, Bernard Cazelles, Luis Aguiar-Conraria, and Maria Joana Soares.}
\references{

Aguiar-Conraria L., and Soares M.J., 2011. 
Business cycle synchronization and the Euro: A wavelet analysis. 
Journal of Macroeconomics 33 (3), 477--489.

Aguiar-Conraria L., and Soares M.J., 2011. 
The Continuous Wavelet Transform: A Primer. 
NIPE Working Paper Series 16/2011.

Aguiar-Conraria L., and Soares M.J., 2012. \code{GWPackage}. Available at \url{http://sites.google.com/site/aguiarconraria/joanasoares-wavelets}; accessed September 4, 2013.

Carmona R., Hwang W.-L., and Torresani B., 1998. 
Practical Time Frequency Analysis. Gabor and Wavelet Transforms with an Implementation in S.
Academic Press, San Diego.

Cazelles B., Chavez M., Berteaux, D., Menard F., Vik J.O., Jenouvrier S., and Stenseth N.C., 2008.
Wavelet analysis of ecological time series.
Oecologia 156, 287--304.

Liu P.C., 1994. 
Wavelet spectrum analysis and ocean wind waves. 
In: Foufoula-Georgiou E., and Kumar P., (eds.), Wavelets in Geophysics, Academic Press, San Diego, 151--166.

Liu Y., Liang X.S., and Weisberg R.H., 2007.
Rectification of the Bias in the Wavelet Power Spectrum. 
Journal of Atmospheric and Oceanic Technology 24, 2093--2102.

Tian, H., and Cazelles, B., 2012. \code{WaveletCo}. Available at \url{http://cran.r-project.org/src/contrib/Archive/WaveletCo/}, archived April 2013; accessed July 26, 2013.

Torrence C., and Compo G.P., 1998. 
A practical guide to wavelet analysis. 
Bulletin of the American Meteorological Society 79 (1), 61--78.

Veleda D., Montagne R., and Araujo M., 2012. 
Cross-Wavelet Bias Corrected by Normalizing Scales. 
Journal of Atmospheric and Oceanic Technology 29, 1401--1408.
}

\keyword{ package }
% \seealso{ }
% \examples{ }
