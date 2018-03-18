\name{WaveletComp-package}
\alias{WaveletComp-package}
\alias{WaveletComp}
\docType{package}
\title{Computational Wavelet Analysis}

\description{Wavelet analysis and reconstruction of time series, cross-wavelets and phase difference (with filtering options), 
significance with bootstrap algorithms.}

\details{
\tabular{ll}{
Package: \tab WaveletComp\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2018-03-18\cr
License: \tab GPL-2\cr
URL: \tab Guide booklet at \cr
     \tab \url{http://www.hs-stat.com/projects/WaveletComp/WaveletComp_guided_tour.pdf} \cr
}

Periodic phenomena of a single time series can be analyzed with function \code{analyze.wavelet}.
Results of the analysis (a time/period image of the wavelet power spectrum, plots of the average power, 
and phase plots for selected periods and a time/period image of phases) can be accessed through various plot functions
(\code{wt.image}, \code{wt.avg}, \code{wt.sel.phases}, \code{wt.phase.image}). 
Function \code{reconstruct} returns the reconstructed time series where reconstruction is according to constraints 
on significance, period specification, and cone of influence.

The cross-wavelet spectrum and coherency spectrum of two time series can be analyzed with function \code{analyze.coherency}.
Results (a time/period image of cross-wavelet power or coherency, plots of averages, plots of phases and phase differences 
for selected periods and the time/period image of phase differences) can be accessed through corresponding functions
(\code{wc.image}, \code{wc.avg}, \code{wc.sel.phases}, \code{wc.phasediff.image}).

Detrending of the time series at hand is offered as an option. Wavelet transformations are computed using the Morlet wavelet. 
Smoothing filters are provided in the case of cross-wavelet transformation to compute wavelet coherency.

Significance is assessed with simulation algorithms, a variety of alternative hypotheses to test is available, for which 
surrogate time series are provided: white noise, shuffling the given time series, time series with a similar spectrum, AR, and ARIMA.

Names and parts of the layout of some routines were inspired by similar functions developed by Huidong Tian 
and Bernard Cazelles (archived R package \code{WaveletCo}). The basic concept of the simulation algorithm 
and of ridge determination build on ideas developed by these authors. The major part of the code for the computation
of the cone of influence and the code for Fourier-randomized surrogate time series has been adopted from Huidong Tian. 
The implementation of a choice of filtering windows for the computation of the wavelet coherence was inspired by 
Luis Aguiar-Conraria and Maria Joana Soares (\code{GWPackage}).

Cross-wavelet and coherence computation, the simulation algorithm, and ridge determination build heavily on the use of matrices
in order to minimize computation time in R. \cr

What is new in \code{WaveletComp} version 1.1?

Tools for displaying and analyzing periodic phenomena across time have been extended. The main innovations are:
\itemize{
\item All functions of family \code{wt.<>} (showing results concerning a single time series) can now also be applied 
to extract univariate outcomes from cross-wavelet and coherence analysis (objects of class \code{"analyze.coherency"}).
\item It is possible to control the color gradation of time-period spectrum plots, and accentuate the contrast, by raising 
the wavelet power values to any (positive) exponent before plotting.
\item Setting a maximum level for the color bar facilitates the visual comparison of time-period spectrum plots. 
Maximum and minimum plot levels are options for plots of averages too.
\item The time and period axes are now easier to individualize by specifying tick marks and labels. 
Coordinates on the time axis can be conveniently addressed via an index or a \code{"POSIXct"} object.
\item Graphical parameters of global coverage (\code{cex.axis}, \code{font.axis}, \code{cex.lab}, \code{font.lab}, 
\code{mgp} etc., see \code{par}) as well as parameters of local coverage (within axis specification options) help fine-tune plots.
\item Two more real-world data sets have been included in WaveletComp, namely:
\itemize{
\item  Data set \code{"weather.radiation.Mannheim"}, containing daily weather and ambient radiation readings from Mannheim (Germany).
\item Data set \code{"USelection2016.Instagram"}, containing hourly numbers of candidate-related media uploads to Instagram right before the 2016 US presidential election.
}
}
}
\author{Angi Roesch and Harald Schmidbauer; credits are also due to Huidong Tian, Bernard Cazelles, Luis Aguiar-Conraria, and Maria Joana Soares.}
\references{

Aguiar-Conraria L., and Soares M.J., 2011. 
Business cycle synchronization and the Euro: A wavelet analysis. 
Journal of Macroeconomics 33 (3), 477--489.

Aguiar-Conraria L., and Soares M.J., 2011. 
The Continuous Wavelet Transform: A Primer. 
NIPE Working Paper Series 16/2011.

Aguiar-Conraria L., and Soares M.J., 2012. \code{GWPackage}. 
Available at \url{https://sites.google.com/site/aguiarconraria/joanasoares-wavelets}; accessed September 4, 2013.

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

Schmidbauer H., Roesch A., Stieler F., 2018.
The 2016 US presidential election and media on Instagram: Who was in the lead?
Computers in Human Behavior 81, 148--160. 
\doi{10.1016/j.chb.2017.11.021} 

Tian, H., and Cazelles, B., 2012. \code{WaveletCo}. 
Available at \url{https://cran.r-project.org/src/contrib/Archive/WaveletCo/}, archived April 2013; accessed July 26, 2013.

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
