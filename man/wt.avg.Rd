\name{wt.avg}
\alias{wt.avg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of wavelet power averages across time of a single time series
}
\description{
This function plots wavelet power averages across time of a single time
series, which are provided by an object of class \code{analyze.wavelet}, or alternatively of class \code{analyze.coherency}. 
(In the latter case, the series number or name can be specified.) 
The vertical axis shows the Fourier periods. The horizontal axis shows the averages.

There is an option to label periods according to significance of
averages (if p-values are provided by \code{analyze.wavelet}) at
given levels of significance. Labels are point symbols along the
line of averages which can be assigned individually.

The idea to show significance levels by colors of plotting characters and its implementation has been adopted from Huidong Tian and Bernard Cazelles 
(archived R package \code{WaveletCo}).
}
\usage{
wt.avg(WT, my.series = 1, 
  show.siglvl = T, siglvl = c(0.05, 0.1), sigcol = c("red", "blue"), sigpch = 20, 
  label.avg.axis = T, averagelab = NULL, 
  label.period.axis = T, periodlab = NULL, 
  show.legend = T, legend.coords = "topright", 
  main = NULL, lwd = 0.5,
  verbose = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{WT}{an object of class \code{analyze.wavelet} or \code{analyze.coherency}}
  \item{my.series}{In case class(WT) = \code{analyze.coherency}: number (1 or 2) or name of the series to be analyzed. Default: 1. }
  
  \item{show.siglvl}{Label periods according to significance of averages? (Effective only if p-values are provided by \code{analyze.coherency}.) 
  
                     Default: TRUE.
  }
  \item{siglvl}{a vector of significance levels (of any length and order). Default: c(0.05, 0.1).}
  \item{sigcol}{a vector of colors (should be of same length as and correspond to \code{siglvl}, otherwise 1:length(siglvl)).
  
                Default: c("red","blue").
  }
  \item{sigpch}{a vector of plotting characters. (It should be of same length as and correspond to \code{siglvl} to produce different labels, otherwise the default setting is used. A single input value affects all labels.) Default: 20.}
  \item{label.avg.axis}{Label the axis of averages? Logical. Default: TRUE.}
  \item{averagelab}{Label for the axis of averages.
  
                    Default: "average wavelet power".
  }
  \item{label.period.axis}{Label the (Fourier) period axis? Logical. Default: TRUE.}
  \item{periodlab}{(Fourier) period axis label. Default: "period".} 
  \item{show.legend}{Include legend of significance levels into the plot? Logical. Default: TRUE.}
  \item{legend.coords}{coordinates to position the legend (as in function \code{legend}). Default: "topright".}
  \item{main}{an overall title for the plot. Default: NULL.}
  \item{lwd}{line width. Default: 0.5.}
  \item{verbose}{Print verbose output on the screen? Logical. Default: FALSE.}
}
% \details{ }
% \value{ }
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
\author{Angi Roesch and Harald Schmidbauer; credits are also due to Huidong Tian and Bernard Cazelles}
% \note{ }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ 
\code{\link{analyze.wavelet}}, \code{\link{wt.image}}, \code{\link{wt.sel.phases}}, \code{\link{wt.phase.image}}, \code{\link{reconstruct}} 
}
\examples{
\dontrun{
## The following example is adopted from Liu et al, 2007:

series.length = 6*128*24
x1 = periodic.series(start.period = 1*24, length = series.length)
x2 = periodic.series(start.period = 8*24, length = series.length)
x3 = periodic.series(start.period = 32*24, length = series.length)
x4 = periodic.series(start.period = 128*24, length = series.length)
x = x1 + x2 + x3 + x4
plot(ts(x, start=0, frequency=24), type="l", 
        xlab="time (days)", ylab="hourly data", 
        main="a series of hourly data with periods of 1, 8, 32, and 128 days")
        
my.data = data.frame(x=x)

my.w = analyze.wavelet(my.data, "x", loess.span=0, dt=1/24, dj=1/20, 
                       lowerPeriod=1/4, make.pval=T, n.sim=10)

## Plot of wavelet power spectrum (with equidistant color breakpoints):  
wt.image(my.w, color.key="i", 
         legend.params=list(lab="wavelet power levels (equidistant levels)"))
## Plot of average wavelet power:
wt.avg(my.w, siglvl=0.05, sigcol="red")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }