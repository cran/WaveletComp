\name{wt.image}
\alias{wt.image}

\title{
Image plot of the wavelet power spectrum of a single time series
}
\description{
This function plots the wavelet power spectrum of a single time series, which is provided by an object of class \code{analyze.wavelet}, or alternatively of class
\code{analyze.coherency}. (In the latter case,
the series number or name can be specified.)
The vertical axis shows the Fourier
periods. The horizontal axis shows time step counts, but can
be easily transformed into a calendar axis if dates are provided in
either rownames or as a variable named date in the data frame at
hand. Both axes can be relabeled.

The color levels can be defined according to quantiles of power or
according to equidistant breakpoints (covering the interval from 0
to maximum power), with the number of levels as a further
parameter. In addition, there is an option to adopt an individual
color palette.

Further plot design options concern: plot of the cone of
influence, plot of wavelet power contour lines at a specified
level of significance, plot of power ridges.

Finally, there is an option to insert and format a color legend (a
right-hand vertical color bar) and to set the plot title.  For
further processing of the plot, graphical parameters of plot
regions are provided as output.

The name and parts of the layout were inspired by a similar function developed by Huidong Tian and Bernard Cazelles 
(archived R package \code{WaveletCo}).
}
\usage{
wt.image(WT, my.series = 1, 
     plot.coi = T, 
     plot.contour = T, siglvl = 0.1, col.contour = "white",
     plot.ridge = T, lvl = 0, col.ridge = "black", 
     color.key = "quantile", 
     n.levels = 100, color.palette = "rainbow(n.levels, start=0, end=.7)",
     useRaster = T, max.contour.segments = 250000,
     plot.legend = T,
     legend.params = list(width = 1.2, shrink = 0.9, mar = 5.1, 
                          n.ticks = 6, label.digits = 1, label.format = "f", 
                          lab = NULL, lab.line = 2.5),
     label.time.axis = T, show.date = F, date.format = NULL, timelab = NULL, 
     label.period.axis = T, periodlab = NULL, 
     main = NULL,
     lwd = 2,
     graphics.reset = T,
     verbose = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{WT}{an object of class \code{analyze.wavelet} or \code{analyze.coherency}}
  \item{my.series}{In case class(WT) = \code{analyze.coherency}: number (1 or 2) or name of the series to be analyzed. Default: 1. }
  \item{plot.coi}{Plot cone of influence? Logical. Default: TRUE}
  \item{plot.contour}{Plot contour lines to border the area of wavelet power significance at level \code{siglvl}? Logical. Default: TRUE.}
  \item{siglvl}{level of wavelet power significance applied to the plot of contour lines. Default: 0.1.}
  \item{col.contour}{Color of contour lines. Default: "white".}
  \item{plot.ridge}{Plot the wavelet power ridge? Logical. Default: TRUE.}
  \item{lvl}{minimum level of wavelet power for the ridge to be plotted. Default: 0.}
  \item{col.ridge}{color of the power ridge. Default: "black".}
  \item{color.key}{How to assign colors to power levels? Two options:
  
                   \tabular{rlll}{
                            \tab  "interval" or "i" \tab : \tab equidistant breakpoints \cr
                            \tab                    \tab   \tab (from 0 through maximum value) \cr
                            \tab  "quantile" or "q" \tab : \tab quantiles 
                           } 
                       
                   Default: "quantile"
  
                                    
  }
  \item{n.levels}{Number of color levels. Default: 100.}
  \item{color.palette}{Definition of color levels. (It will be assigned to levels in reverse order!) Default: "rainbow(n.levels, start=0, end=.7)".}
  \item{useRaster}{Use a bitmap raster instead of polygons to plot the wavelet power image? Logical. Default: TRUE.}
  \item{max.contour.segments}{limit on the number of segments in a single contour line, positive integer. Default: 250000 (\code{options(...)} default settings: 25000)}
  \item{plot.legend}{Plot color legend (a vertical bar of colors and breakpoints)? Logical. Default: TRUE}
  \item{legend.params}{a list of parameters for the plot of color legend, parameter values can be set selectively (style in parts adopted from image.plot in the R package "fields" by Douglas Nychka):
  
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
                                \tab               \tab   \tab Default: 1. \cr
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
  \item{show.date}{Show calendar dates? (Effective only if dates are available as rownames or as variable \code{date} in the data frame analyzed using \code{analyze.wavelet}.) Logical. Default: FALSE.}
  \item{date.format}{the format of date given as a character string, e.g. "\%Y-\%m-\%d", or equivalently "\%F"; see \code{strptime} for a list of implemented date conversion specifications. (If not specified, as.Date will be applied.) Default: NULL.}
  \item{timelab}{Time axis label. Default: "time".}
  \item{label.period.axis}{Label the (Fourier) period axis? Logical. Default: TRUE.}
  \item{periodlab}{(Fourier) period axis label. Default: "period".}
  \item{main}{an overall title for the plot. Default: NULL.}
  \item{lwd}{line width of contour lines and ridge. Default: 2.}
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
\code{\link{analyze.wavelet}}, \code{\link{wt.avg}}, \code{\link{wt.sel.phases}}, \code{\link{wt.phase.image}}, \code{\link{reconstruct}} 
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
     
my.date =  seq(as.POSIXct("2014-10-14 00:00:00","\%F \%T"), by="hour", 
               length.out=series.length)     
my.data = data.frame(date=my.date, x=x)

my.w = analyze.wavelet(my.data, "x", 
                       loess.span=0, 
                       dt=1/24, dj=1/20, 
                       lowerPeriod=1/4, 
                       make.pval=T, n.sim=10)

## Plot of wavelet power spectrum with breakpoints referring to quantiles:  
wt.image(my.w, 
   legend.params=list(lab="wavelet power levels (quantiles)", lab.line=3.5, 
                      label.digits=2))

## The same plot, but with equidistant color breakpoints: 
wt.image(my.w, color.key="i", 
         legend.params=list(lab="wavelet power levels (equidistant levels)"))

## The plot with calendar axis:
wt.image(my.w, color.key="i", 
         legend.params=list(lab="wavelet power levels (equidistant levels)"), 
         show.date=T, date.format="\%F \%T", timelab="")

## The same plot, but now with individualized calendar axis:
my.plot = wt.image(my.w, color.key="i", 
                   legend.params=list(lab="wavelet power levels (equidistant levels)"), 
                   label.time.axis=F)
## recover plot region:
par(new=T, plt=my.plot$image.plt)
## empty plot
plot(my.date, rep(1,series.length), type="n", 
     xaxs = "i", yaxs="i", xaxt="n", yaxt="n", 
     xlab="", ylab="")
## individualized calendar axis:
axis.POSIXct(1, at=
 seq(as.POSIXct("2014-11-01 00:00:00", "\%F \%T"), my.date[length(my.date)], by="month"), 
 format="\%b \%Y", las=2)
## return to default plot region:
par(my.plot$op)

## A plot with different colors:
wt.image(my.w, 
   legend.params=list(lab="wavelet power levels (quantiles)", lab.line=3.5, 
                      label.digits=2), 
   color.palette="gray((1:n.levels)/n.levels)", col.ridge="yellow")

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }