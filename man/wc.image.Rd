\name{wc.image}
\alias{wc.image}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Image plot of the cross-wavelet power spectrum and wavelet coherence
spectrum of two time series
}
\description{
This function plots the cross-wavelet power image, or alternatively the wavelet
coherence image, of two time series, which are provided by an object of class
\code{analyze.coherency}.  The vertical axis shows the Fourier
periods. The horizontal axis shows time step counts, but it can
be easily transformed into a calendar axis if dates are provided in
either rownames or a variable named date in the data frame at
hand. Both axes can be relabeled.

An option is given to raise values by any exponent before
plotting.

The color levels can be defined according to quantiles of values
or according to equidistant breakpoints (covering the interval
from 0 to maximum value), with the number of levels as a further
parameter. In addition, there is an option to adopt an individual
color palette.

Further plot design options concern: plot of the cone of
influence, plot of contour lines to border areas of significance, plot of 
the ridge, and plot of arrows (optional: "smoothed" arrows computed 
from smoothing filters as defined in \code{analyze.coherency}) to reflect
phase differences.

For that matter, the significance level of contour lines can be
defined separately. The plot of the ridge can be restricted to a high-level region 
("high" according to a given level of plotted values). 
In particular, the area to be filled with arrows can be determined in several ways: 
to reflect significance (at a given level) with respect to cross-wavelet power, wavelet
coherency, or individual wavelet power, and/or to flag a high-value region. 
Furthermore, there is an option to clear out the area where the p-values of
cross-wavelet power (coherence, respectively) exceed a given
level.

Finally, there is an option to format and insert a color legend (a
right-hand vertical color bar) and to set the plot title.  For
further processing of the plot, graphical parameters of plot
regions are provided as output.

The name and parts of the layout were inspired by a similar function developed by Huidong Tian and Bernard Cazelles 
(archived R package \code{WaveletCo}). 
The code for the arrow design to reflect phase differences 
has been adopted from Huidong Tian.
}
\usage{
wc.image(WC, 
     which.image = "wp", exponent = 1, 
     plot.coi = T, 
     plot.contour = T, siglvl.contour = 0.1, col.contour = "white",
     plot.ridge = F, lvl = 0, col.ridge = "black",
     plot.arrow = T, use.sAngle = F, 
     p = 1, 
     which.arrow.sig = which.image, siglvl.arrow = 0.05, col.arrow = "black",
     clear.area = F, 
     which.area.sig = which.image, siglvl.area = 0.2, 
     color.key = "quantile", 
     n.levels = 100, color.palette = "rainbow(n.levels, start=0, end=.7)", 
     useRaster = T, max.contour.segments = 250000,
     plot.legend = T,
     legend.params = list(width=1.2, shrink=0.9, mar=5.1, 
                          n.ticks=6, label.digits=1, label.format="f", 
                          lab=NULL, lab.line=2.5),
     label.time.axis = T, show.date = F, date.format = NULL, timelab = NULL, 
     label.period.axis = T, periodlab = NULL, 
     main = NULL,
     lwd = 2,
     graphics.reset = T,
     verbose = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{WC}{an object of class \code{analyze.coherency}}
  \item{which.image}{Which image is to be plotted?
                     
                       \tabular{rlll}{
                                \tab "wp" \tab : \tab cross-wavelet power \cr
                                \tab "wc" \tab : \tab wavelet coherence \cr
                                }
                                
                        Default: "wp"
 
  }
  \item{exponent}{Exponent of values to be plotted. Default: 1.}
  \item{plot.coi}{Plot cone of influence? Logical. Default: TRUE}
  \item{plot.contour}{Plot contour lines to border the area of cross-wavelet power, respectively wavelet coherence significance at level \code{siglvl.contour}? Logical. Default: TRUE.}
  \item{siglvl.contour}{level of cross-wavelet power, respectively wavelet coherence significance applied to the plot of contour lines. Default: 0.1.}
  \item{col.contour}{Color of contour lines. Default: "white".}
  \item{plot.ridge}{Plot the cross-wavelet, resp. wavelet coherence power ridge? Logical. Default: FALSE.}
  \item{lvl}{minimum level of cross-wavelet power (or wavelet coherence) for the ridge to be plotted, or alternatively, within the area of arrows (if p=0 or 2). Default: 0.}
  \item{col.ridge}{color of the cross-wavelet power, resp. wavelet coherence ridge. Default: "black".}
  \item{plot.arrow}{Plot arrows depicting the phase difference? Logical. Default: TRUE.}
  \item{use.sAngle}{Use smoothed version of phase difference? Logical. Default: FALSE.}
  \item{p}{Which area should be filled with arrows displaying phase differences?\cr (Only effective if \code{plot.arrow}=TRUE.)
  
           \tabular{rlll}{
                       \tab p=0 \tab : \tab area with high values of \code{which.lvl} only \cr
                       \tab     \tab   \tab (cf. \code{lvl}) \cr
                       \tab p=1 \tab : \tab area of significance of \code{which.sig} only \cr
                       \tab     \tab   \tab (cf. \code{siglvl}) \cr
                       \tab p=2 \tab : \tab area with both high values and significance \cr
                       \tab     \tab   \tab (combining p=0 and p=2)
                       }
    
           Default: 1
  }  
  \item{which.arrow.sig}{Which spectrum and corresponding p-values should be used to restrict the area of arrows according to significance?
                           
                       \tabular{rlll}{
                                \tab "wp" \tab : \tab cross-wavelet power \cr
                                \tab "wc" \tab : \tab wavelet coherence \cr
                                \tab "wt" \tab : \tab individual wavelet power
                                }
                                
                        Default: \code{which.image}

                        }
  \item{siglvl.arrow}{level of significance referring to \code{which.arrow.sig} (if plot.arrow=TRUE and p=1 or 2).

                      Default: 0.05
  }
  \item{col.arrow}{color for the plot of arrows. Default: "black".}
  \item{clear.area}{Clear out an area where p-values are above a certain level? Logical.
  
                    (p-values will refer to the spectrum defined by which.area.sig and significance level siglvl.area.) 
                    
                    Default: FALSE
  }
  \item{which.area.sig}{Which power spectrum and corresponding p-values should be used to clear the outer area? (if clear.area=TRUE) 
  
                        \tabular{rlll}{
                                \tab "wp" \tab : \tab cross-wavelet power \cr
                                \tab "wc" \tab : \tab wavelet coherence \cr
                                \tab "wt" \tab : \tab individual wavelet power
                                }
                                
                        Default: \code{which.image}

                        }
  \item{siglvl.area}{level of significance referring to which.area.sig (if clear.area=TRUE) 

                     Default: 0.2
  }
  \item{color.key}{How to assign colors to power and coherence levels? Two options:
  
                   \tabular{rlll}{
                            \tab  "interval" or "i" \tab : \tab equidistant breakpoints \cr
                            \tab                    \tab   \tab (from 0 through maximum value) \cr
                            \tab  "quantile" or "q" \tab : \tab quantiles 
                           } 
                       
                   Default: "quantile"
  }
  \item{n.levels}{Number of color levels. Default: 100.}
  \item{color.palette}{Definition of color levels. (It will be assigned to levels in reverted order!) Default: "rainbow(n.levels, start=0, end=.7)".}
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
  \item{show.date}{Show calendar dates? (Effective only if dates are available as rownames or by variable \code{date} in the data frame which has been analyzed using  \code{analyze.coherency}.) Logical. Default: FALSE.}
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
\author{Angi Roesch and Harald Schmidbauer; credits are also due to Huidong Tian, and Bernard Cazelles.}
% \note{ }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ 
\code{\link{analyze.coherency}}, \code{\link{wc.avg}}, \code{\link{wc.sel.phases}}, \code{\link{wc.phasediff.image}}
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

my.date = seq(as.POSIXct("2014-10-14 00:00:00","\%F \%T"), by="hour", 
              length.out=series.length)  
my.data = data.frame(date=my.date, x=x, y=y)

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

## plot of cross-wavelet power, with color breakpoints according to quantiles:
wc.image(my.wc, timelab="time (days)", periodlab="period (days)", 
         main="cross-wavelet power",
         legend.params=list(lab="cross-wavelet power levels (quantiles)"))

## The same plot, but with equidistant color breakpoints: 
wc.image(my.wc, color.key="i", timelab="time (days)", periodlab="period (days)", 
         main="cross-wavelet power",
         legend.params=list(lab="cross-wavelet power levels (equidistant levels)"))

## The same plot, but adopting a palette of gray colors:
wc.image(my.wc, color.key="i", timelab="time (days)", periodlab="period (days)", 
         main="cross-wavelet power", 
         legend.params=list(lab="cross-wavelet power levels (equidistant levels)"),
         color.palette="gray( (1:n.levels)/n.levels )", plot.arrow=F)
         
## The same plot, but with yellow arrows and calendar axis:
wc.image(my.wc, color.key="i", timelab="", periodlab="period (days)", 
         main="cross-wavelet power", 
         legend.params=list(lab="cross-wavelet power levels (equidistant levels)"),
         color.palette="gray( (1:n.levels)/n.levels )", 
         col.arrow="yellow", 
         show.date=T)
         
## With additional ridge:
wc.image(my.wc, color.key="i", timelab="", periodlab="period (days)", 
         main="cross-wavelet power", 
         legend.params=list(lab="cross-wavelet power levels (equidistant levels)"),
         color.palette="gray( (1:n.levels)/n.levels )", 
         col.arrow="yellow", 
         show.date=T,
         plot.ridge=T, col.ridge="red")         
         
         
## The same plot, but with yellow arrows and individualized calendar axis:
my.plot = wc.image(my.wc, color.key="i", timelab="", periodlab="period (days)", 
             main="cross-wavelet power", 
             legend.params=list(lab="cross-wavelet power levels (equidistant levels)"),
             color.palette="gray( (1:n.levels)/n.levels )", 
             col.arrow="yellow", 
             label.time.axis =F)                   
## recover plot region:
par(new=T, plt=my.plot$image.plt)
## empty plot
plot(my.date, rep(1,series.length), type="n",
     xaxs = "i", yaxs ="i", xaxt="n", yaxt="n", 
     xlab="", ylab="")
## individualized calendar axis:
axis.POSIXct(1, at=
 seq(as.POSIXct("2014-11-01 00:00:00", "\%F \%T"), my.date[length(my.date)], by="month"), 
 format="\%b \%Y", las=2)
## return to default plot region:
par(my.plot$op)         

## plot of wavelet coherence, with color breakpoints according to quantiles:
wc.image(my.wc, which.image="wc", 
  timelab="time (days)", periodlab="period (days)", 
  main="wavelet coherence",
  legend.params=list(lab="wavelet coherence levels (quantiles)", lab.line=3.5, 
                     label.digits=3))
## plot of wavelet coherence, but with equidistant color breakpoints:
wc.image(my.wc, which.image="wc", color.key="i", 
         timelab="time (days)", periodlab="period (days)", 
         main="wavelet coherence",
         legend.params=list(lab="wavelet coherence levels (equidistant levels)"))         
         
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }