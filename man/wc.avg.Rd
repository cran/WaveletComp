\name{wc.avg}
\alias{wc.avg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot cross-wavelet power averages and wavelet coherence averages across time of two time series
}
\description{
This function plots cross-wavelet power averages across time, or alternatively 
wavelet coherence averages, of two time series, 
which are provided by an object of class \code{analyze.coherency}. 
The vertical axis shows the Fourier periods. The horizontal axis shows the averages.

There is an option to label periods according to significance of
averages (if p-values are provided by \code{analyze.coherency}) at
given levels of significance. Labels are point symbols along the
line of averages which can be assigned individually.

The idea to show significance levels by colors of plotting characters and its implementation has been adopted from Huidong Tian and Bernard Cazelles 
(archived R package \code{WaveletCo}).
}
\usage{
wc.avg(WC, which.avg = "wp", show.siglvl = T, siglvl = c(0.05, 0.1), 
   sigcol = c("red", "blue"), sigpch = 20, 
   label.avg.axis = T, averagelab = NULL, 
   label.period.axis = T, periodlab = NULL, 
   show.legend = T, legend.coords = "topright", 
   main = NULL, lwd = 0.5,
   verbose = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{WC}{an object of class \code{analyze.coherency}.}
  \item{which.avg}{Which averages should be plotted?
  
                   \tabular{rlll}{
                            \tab "wp" \tab : \tab cross-wavelet power \cr
                            \tab "wc" \tab : \tab wavelet coherence                            
                           }
                                
                   Default: "wp"                
  }
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
  
                    Default: "average cross-wavelet power" (if which.avg="wp"), "average coherence" (if which.avg="wc").
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
Business cycle synchronization and the Euro: A wavelet analysis. 
Journal of Macroeconomics 33 (3), 477--489.

Aguiar-Conraria L., and Soares M.J., 2011. 
The Continuous Wavelet Transform: A Primer. 
NIPE Working Paper Series 16/2011.

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
\author{Angi Roesch and Harald Schmidbauer; credits are also due to Huidong Tian and Bernard Cazelles.}
% \note{ }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{analyze.coherency}}, \code{\link{wc.image}}, \code{\link{wc.sel.phases}}, \code{\link{wc.phasediff.image}}
}
\examples{
\dontrun{
## The following example is adopted from Veleda et al, 2012:

add.noise=TRUE

series.length = 3*128*24
x1 = periodic.series(start.period = 1*24, length = series.length)
x2 = periodic.series(start.period = 2*24, length = series.length)
x3 = periodic.series(start.period = 4*24, length = series.length)
x4 = periodic.series(start.period = 8*24, length = series.length)
x5 = periodic.series(start.period = 16*24, length = series.length)
x6 = periodic.series(start.period = 32*24, length = series.length)
x7 = periodic.series(start.period = 64*24, length = series.length)
x8 = periodic.series(start.period = 128*24, length = series.length)

x = x1 + x2 + x3 + x4 + 3*x5 + x6 + x7 + x8 
y = x1 + x2 + x3 + x4 + 3*x5 + x6 + 3*x7 + x8

if (add.noise == TRUE){
    x = x + rnorm(length(x))
    y = y + rnorm(length(y))
}

my.data = data.frame(x=x, y=y)

ts.plot(ts(my.data$x, start=0, frequency=24), 
     ts(my.data$y, start=0, frequency=24),
     type="l", col=1:2, 
     xlab="time (days)", ylab="hourly data", 
     main="a series of hourly data with periods of 1, 2, 4, 8, 16, 32, 64, and 128 days",
     sub="(different amplitudes at periods 16 and 64)")
legend("topright", legend=c("x","y"), col=1:2, lty=1)

## computation of cross-wavelet power and wavelet coherence:
my.wc = analyze.coherency(my.data, c("x","y"), loess.span=0, 
                           dt=1/24, dj=1/20, 
                           window.size.t=1, window.size.s=1/2, 
                           lowerPeriod=1/4, 
                           make.pval=T, n.sim=10)

## plot of cross-wavelet power (with color breakpoints according to quantiles):
wc.image(my.wc, timelab="time (days)", periodlab="period (days)", 
         main="cross-wavelet power")
## plot of average cross-wavelet power:
wc.avg(my.wc, siglvl=0.05, sigcol="red")

## plot of wavelet coherence (with color breakpoints according to quantiles):
wc.image(my.wc, which.image="wc", timelab="time (days)", periodlab="period (days)", 
         main="wavelet coherence",
         legend.params=list(label.digits=3))
## plot of average wavelet coherence:
wc.avg(my.wc, which.avg="wc", siglvl=0.05, sigcol="red", legend.coords="topleft")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }