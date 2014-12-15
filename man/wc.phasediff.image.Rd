\name{wc.phasediff.image}
\alias{wc.phasediff.image}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Image plot of phase differences of periodic components for two time
series
}
\description{
This function plots the phase difference image of two time series, which is provided by an object of class \code{analyze.coherency}.  
The vertical axis shows the Fourier periods. The horizontal axis shows time step counts,
but can be easily transformed into a calendar axis if
dates are provided in either rownames or a variable named date in
the data frame at hand. Both axes can be relabeled.

The color levels are defined according to equidistant breakpoints
(covering the interval from -pi to +pi), with the number of levels
as a further parameter. In addition, there is an option to adopt
an individual color palette.

If the default palette is retained, colors indicate the following. Green:
phase differences close to zero, which means that the two time
series are in phase at the respective period. Yellowgreen: in
phase, series 1 leading. Turquoise: in phase, series 2 leading.
Red: phase differences are close to +pi, out of phase, series
2 leading. Blue: phase differences are close to -pi, out of
phase, series 1 leading.

Further plot design options concern: plot of the cone of
influence, plot of contour lines to border areas of significance
with respect to cross-wavelet power or wavelet coherence at a
given significance level.

Finally, there is an option to insert and format a color legend (a
right-hand vertical color bar) and to set the plot title.  For
further processing of the plot, graphical parameters of plot
regions are provided as output.
}
\usage{
wc.phasediff.image(WC, use.sAngle = F, 
   plot.coi = T, 
   plot.contour = T, which.contour = "wp", siglvl = 0.1, col.contour = "white",
   n.levels = 100, color.palette = "rainbow(n.levels, start=0, end=.7)", 
   useRaster = T, max.contour.segments = 250000,
   plot.legend = T,
   legend.params = list(width=1.2, shrink=0.9, mar=5.1, 
                        n.ticks=6, label.digits=2, label.format="f", 
                        lab=NULL, lab.line=2.5),
   label.time.axis = T, show.date = F, date.format = NULL, timelab = NULL, 
   label.period.axis = T, periodlab = NULL, 
   main = NULL,
   graphics.reset = T,
   verbose = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{WC}{an object of class \code{analyze.coherency}}
  \item{use.sAngle}{Use smoothed version of phase difference? Logical. Default: FALSE.}
  \item{plot.coi}{Plot cone of influence? Logical. Default: TRUE}
  \item{plot.contour}{Plot contour lines to border the area of cross-wavelet power, respectively wavelet coherence significance at level siglvl.contour? Logical. Default: TRUE.}
  \item{which.contour}{Contour lines of which spectrum should be plotted?
  
                       \tabular{rlll}{
                                \tab "wp" \tab : \tab cross-wavelet power \cr
                                 \tab "wc" \tab : \tab wavelet coherence 
                                }
                                
                        Default: "wp"
                       }
  \item{siglvl}{level of cross-wavelet power, respectively wavelet coherence significance applied to the plot of contour lines. Default: 0.1.}
  \item{col.contour}{Color of contour lines. Default: "white".}
  \item{n.levels}{Number of color levels. Default: 100.}
  \item{color.palette}{Definition of color levels. (It will be assigned to levels in reverse order!) Default: "rainbow(n.levels, start=0, end=.7)".}
  \item{useRaster}{Use a bitmap raster instead of polygons to plot the wavelet power image? Logical. Default: TRUE.}
  \item{max.contour.segments}{limit on the number of segments in a single contour line, positive integer. Default: 250000 (\code{options(...)} default settings: 25000)}
  \item{plot.legend}{Plot color legend (a vertical bar of colors and breakpoints)? Logical. Default: TRUE}
  \item{legend.params}{a list of parameters for the plot of color legend, parameter values can be set selectively (style in parts adopted from image.plot in the R package fields by Douglas Nychka):
  
                       \tabular{rlll}{
                                \tab width         \tab : \tab width of legend bar. \cr 
                                \tab               \tab   \tab Default: 1.2. \cr
                                \tab shrink        \tab : \tab a vertical shrinkage factor. \cr
                                \tab               \tab   \tab Default: 0.9. \cr
                                \tab mar           \tab : \tab right margin of legend bar. \cr
                                \tab               \tab   \tab Default: 5.1. \cr
                                \tab n.ticks       \tab : \tab number of ticks for labels. \cr 
                                \tab               \tab   \tab Default: 6. \cr
                                \tab label.digits  \tab : \tab digits of labels. \cr
                                \tab               \tab   \tab Default: 2. \cr
                                \tab label.format  \tab : \tab format of labels. \cr
                                \tab               \tab   \tab Default: "f". \cr
                                \tab lab           \tab : \tab axis label. \cr
                                \tab               \tab   \tab Default: NULL. \cr
                                \tab lab.line      \tab : \tab line (in user coordinate units) where  \cr
                                \tab               \tab   \tab to put the axis label. \cr
                                \tab               \tab   \tab Default: 2.5.
                                }
  }
  \item{label.time.axis}{Label the time axis? Logical. Default: TRUE.}
  \item{show.date}{Show calendar dates? (Effective only if dates are available as rownames or by variable \code{date} in the data frame which has been analyzed using  \code{analyze.coherency}.) Logical. Default: FALSE.}
  \item{date.format}{the format of date given as a character string, e.g. "\%Y-\%m-\%d", or equivalently "\%F"; see \code{strptime} for a list of implemented date conversion specifications. (If not specified, as.Date will be applied.) Default: NULL.}
  \item{timelab}{Time axis label. Default: "time".}
  \item{label.period.axis}{Label the (Fourier) period axis? Logical. Default: TRUE.}
  \item{periodlab}{(Fourier) period axis label. Default: "period".} 
  \item{main}{an overall title for the plot. Default: NULL.}
  \item{graphics.reset}{Reset graphical parameters? Logical. Default: TRUE}
  \item{verbose}{Print verbose output on the screen? Logical. Default: FALSE.}
}
% \details{ }
\value{A list of class \code{graphical parameters} with the following elements:
   \item{op}{original graphical parameters}
   \item{image.plt}{image plot region}
   \item{legend.plt}{legend plot region} 
}
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

Torrence C., and Compo G.P., 1998. 
A practical guide to wavelet analysis.
Bulletin of the American Meteorological Society 79 (1), 61--78.

Veleda D., Montagne R., and Araujo M., 2012. 
Cross-Wavelet Bias Corrected by Normalizing Scales. 
Journal of Atmospheric and Oceanic Technology 29, 1401--1408.
}
\author{Angi Roesch and Harald Schmidbauer}
% \note{ }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ 
\code{\link{analyze.coherency}}, \code{\link{wc.avg}}, \code{\link{wc.image}}, \code{\link{wc.sel.phases}} 
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
my.wc = analyze.coherency(my.data, c("x","y"), 
                          loess.span=0, 
                          dt=1/24, dj=1/20, 
                          window.size.t=1, window.size.s=1/2, 
                          lowerPeriod=1/4, 
                          make.pval=T, n.sim=10)

## plot of phase differences (with contour lines referring to cross-wavelet power)
wc.phasediff.image(my.wc, which.contour="wp", timelab="time (days)", 
                   main="image of phase differences")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }