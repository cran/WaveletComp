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
\code{"analyze.coherency"}.  
The vertical axis shows the Fourier periods. The horizontal axis shows time step counts, but can
be easily transformed into a calendar axis if dates are provided in either row names or a variable 
named \code{"date"} in the data frame at hand. Both axes can be relabeled.
In particular, an option is given to individualize the period and/or time axis 
by specifying tick marks and labels.

An option is given to raise cross-wavelet power (or wavelet coherence) values 
to any (positive) exponent before plotting in order to accentuate the contrast of the image.

The color levels can be defined according to quantiles of values
or according to equidistant breakpoints (covering the interval
from 0 to maximum level), with the number of levels as a further
parameter. A user-defined maximum level can be applied to cross-wavelet power images. 
In addition, there is an option to adopt an individual color palette.

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
coherence, or individual wavelet power, and/or to flag a high-value region. 
Furthermore, there is an option to clear out the area where the p-values of
cross-wavelet power (coherence, respectively) exceed a given
level.

Finally, there is an option to format and insert a color legend (a
right-hand vertical color bar) and to set the plot title.  For
further processing of the plot, graphical parameters of plot
regions are provided as output.

The name and parts of the layout were inspired by a similar function developed by 
Huidong Tian and Bernard Cazelles (archived R package \code{WaveletCo}). 
The code for the arrow design to reflect phase differences 
has been adopted from Huidong Tian.
}
\usage{
wc.image(WC, 
     which.image = "wp", exponent = 1, 
     plot.coi = TRUE, 
     plot.contour = TRUE, siglvl.contour = 0.1, col.contour = "white",
     plot.ridge = FALSE, lvl = 0, col.ridge = "black",
     plot.arrow = TRUE, use.sAngle = FALSE, 
     p = 1, 
     which.arrow.sig = which.image, 
     siglvl.arrow = 0.05, col.arrow = "black",
     clear.area = FALSE, 
     which.area.sig = which.image, siglvl.area = 0.2, 
     color.key = "quantile", 
     n.levels = 100, 
     color.palette = "rainbow(n.levels, start = 0, end = .7)", 
     maximum.level = NULL,
     useRaster = TRUE, max.contour.segments = 250000,
     plot.legend = TRUE,
     legend.params = list(width=1.2, shrink = 0.9, mar = 5.1, 
                          n.ticks = 6, 
                          label.digits = 1, label.format = "f", 
                          lab = NULL, lab.line = 2.5),
     label.time.axis = TRUE, 
     show.date = FALSE, date.format = NULL, date.tz = NULL,
     timelab = NULL, timetck = 0.02, timetcl = 0.5,
     spec.time.axis = list(at = NULL, labels = TRUE, 
                           las = 1, hadj = NA, padj = NA),
     label.period.axis = TRUE, 
     periodlab = NULL, periodtck = 0.02, periodtcl = 0.5,
     spec.period.axis = list(at = NULL, labels = TRUE, 
                             las = 1, hadj = NA, padj = NA),
     main = NULL,
     lwd = 2, lwd.axis = 1, 
     graphics.reset = TRUE,
     verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{WC}{an object of class \code{"analyze.coherency"}}
  \item{which.image}{Which image is to be plotted?
                     
                       \tabular{rlll}{
                                \tab \code{"wp"} \tab : \tab cross-wavelet power \cr
                                \tab \code{"wc"} \tab : \tab wavelet coherence \cr
                                }
                                
  Default: \code{"wp"}.
  }
  %%%%%%%%%%%%%%%%%%%%%%%%%
  \item{exponent}{Exponent applied to values before plotting in order to accentuate the contrast of the image; 
  the exponent should be positive.
  
  Default: \code{1}.
  }
  %%%%%%%%%%%%%%%%%%%%%%%%%
  \item{plot.coi}{Plot cone of influence? Logical. Default: \code{TRUE}.}
  %%%%%%%%%%%%%%%%%%%%%%%%%
  \item{plot.contour}{Plot contour lines to border the area of cross-wavelet power (or wavelet coherence, 
  depending on \code{which.image}) significance? Logical.
  
  Default: \code{TRUE}.
  }
  \item{siglvl.contour}{level of cross-wavelet power (or wavelet coherence, depending on \code{which.image}) 
  significance to be applied to the plot of contour lines. 
  
  Default: \code{0.1}.
  }
  \item{col.contour}{color of contour lines. Default: \code{"white"}.}
  %%%%%%%%%%%%%%%%%%%%%%%%%
  \item{plot.ridge}{Plot the cross-wavelet power (or wavelet coherence, depending on \code{which.image}) ridge? Logical.
  
  Default: \code{FALSE}.
  }
  \item{lvl}{minimum level of cross-wavelet power (or wavelet coherence, depending on \code{which.image}) for ridge or arrows to be plotted.\cr
  (Only effective if \code{plot.ridge = TRUE} or, when setting \code{p = 0} or \code{p = 2}, if \code{plot.arrow = TRUE}.)
  
  Default: \code{0}.
  }
  \item{col.ridge}{ridge color. Default: \code{"black"}.}
  %%%%%%%%%%%%%%%%%%%%%%%%%
  \item{plot.arrow}{Plot arrows depicting the phase difference? Logical.
  
  Default: \code{TRUE}.}
  \item{use.sAngle}{Use smoothed version of phase difference? Logical. 
  
  Default: \code{FALSE}.}
  \item{p}{Which area should be filled with arrows displaying phase differences?\cr 
  (Only effective if \code{plot.arrow = TRUE}.)
  \itemize{
  \item[\code{0} :] Select the area w.r.t. power (or coherence, depending on \code{which.image}.\cr 
  The value of \code{lvl} is the threshold.\cr 
  Set \code{lvl = 0} to plot arrows over the entire area.
  \item[\code{1} :] Select the area w.r.t. significance as specified in \code{which.arrow.sig}.\cr 
  The value of \code{siglvl.arrow} is the threshold.
  \item[\code{2} :] Intersection of both areas above. 
  }
  Default: \code{1}.
  }  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  \item{which.arrow.sig}{Which spectrum (and corresponding p-values) should be used to restrict the plot 
  of arrows according to significance?
                           
                       \tabular{rlll}{
                                \tab \code{"wp"} \tab : \tab cross-wavelet power \cr
                                \tab \code{"wc"} \tab : \tab wavelet coherence \cr
                                \tab \code{"wt"} \tab : \tab individual wavelet power
                                }
                                
  Default: \code{which.image}.
  }
  \item{siglvl.arrow}{level of significance for arrows to be plotted. \cr 
  (Only effective if \code{plot.arrow = TRUE} and \code{p = 1} or \code{p = 2}.)

  Default: \code{0.05}.
  }
  \item{col.arrow}{arrow color. Default: \code{"black"}.}
  %%%%%%%%%%%%%%%%%%%%%%%%%
  \item{clear.area}{Clear out an area where p-values are above a certain level? Logical.
  
  (Here, p-values will refer to the spectrum defined by \code{which.area.sig} and significance level \code{siglvl.area}, see below.) 
                    
  Default: \code{FALSE}.
  }
  \item{which.area.sig}{Which power spectrum (and corresponding p-values) should be used to clear the outer area? \cr
  (Only effective if \code{clear.area = TRUE}) 
  
                        \tabular{rlll}{
                                \tab \code{"wp"} \tab : \tab cross-wavelet power \cr
                                \tab \code{"wc"} \tab : \tab wavelet coherence \cr
                                \tab \code{"wt"} \tab : \tab individual wavelet power
                                }
                                
  Default: \code{which.image}.
  }
  \item{siglvl.area}{level of significance for the area to be cleared out. \cr
  (Only effective if \code{clear.area = TRUE}.) 

  Default: \code{0.2}.
  }
  %%%%%%%%%%%%%%%%%%%%%%%%%
  \item{color.key}{How to assign colors to power and coherence levels? Two options:
  
                   \tabular{rlll}{
                            \tab  \code{"interval"} or \code{"i"} \tab : \tab equidistant breakpoints \cr
                            \tab                                  \tab   \tab (from \code{0} through maximum value) \cr
                            \tab  \code{"quantile"} or \code{"q"} \tab : \tab quantiles 
                           } 
                       
  Default: \code{"quantile"}.
  }
  \item{n.levels}{Number of color levels. Default: \code{100}.}
  \item{color.palette}{Definition of color levels. (The color palette will be assigned to levels in reverse order!)
  
  Default: \code{"rainbow(n.levels, start = 0, end = .7)"}.
  }
  %%%%%%%%%%%%%%%%%%%%%%%%%
  \item{maximum.level}{Maximum plot level of cross-wavelet power considered; only effective in case of equidistant breakpoints (\code{color.key} equaling \code{"i"}).  
  
  Default: \code{NULL} (referring to maximum level observed).\cr 
  Wavelet coherence has maximum plot level \code{1} by definition.}
  %%%%%%%%%%%%%%%%%%%%%%%%%
  \item{useRaster}{Use a bitmap raster instead of polygons to plot the image? Logical.
  
  Default: \code{TRUE}.}
  \item{max.contour.segments}{limit on the number of segments in a single contour line, positive integer. 
  
  Default: \code{250000} (\code{options(...)} default settings: \code{25000}).}
  %%%%%%%%%%%%%%%%%%%%%%%%%
  \item{plot.legend}{Plot color legend (a vertical bar of colors and breakpoints)? Logical. 
  
  Default: \code{TRUE}.}
  \item{legend.params}{a list of parameters for the plot of the color legend; parameter values can be set selectively 
  (style in parts adopted from \code{image.plot} in the R package \code{fields} by Douglas Nychka):
  
  \itemize{
  \item[\code{width}:] width of legend bar. \cr 
  Default: \code{1.2}. 
  \item[\code{shrink}:] a vertical shrinkage factor. \cr
  Default: \code{0.9}. 
  \item[\code{mar}:] right margin of legend bar. \cr
  Default: \code{5.1}. 
  \item[\code{n.ticks}:] number of ticks for labels. \cr 
  Default: \code{6}. 
  \item[\code{label.digits}:] digits of labels. \cr
  Default: \code{1}. 
  \item[\code{label.format}:] format of labels. \cr
  Default: \code{"f"}. 
  \item[\code{lab}:] axis label. \cr
  Default: \code{NULL}. 
  \item[\code{lab.line}:] line (in user coordinate units) where to put the axis label. \cr
  Default: \code{2.5}.
          }
  }
  %%%%%%%%%%%%%%%%%%%%%%%%%
  \item{label.time.axis}{Label the time axis? Logical. 
  
  Default: \code{TRUE}.}
  \item{show.date}{Show calendar dates? (Effective only if dates are available as row names or by variable 
  \code{date} in the data frame which has been analyzed.) Logical. 
  
  Default: \code{FALSE}.}
  %%%%%%%%%%%%%%%%%%%%%%%%%
  \item{date.format}{the format of calendar date given as a character string, e.g. \code{"\%Y-\%m-\%d"}, or equivalently \code{"\%F"}; 
  see \code{strptime} for a list of implemented date conversion specifications. Explicit information given here will overturn any specification 
  stored in \code{WC}. If unspecified, date formatting is attempted according to \code{as.Date}. 
  
  Default: \code{NULL}.}
  \item{date.tz}{a character string specifying the time zone of calendar date; see \code{strptime}. Explicit information given here will overturn 
  any specification stored in \code{WC}. If unspecified, \code{""} (the local time zone) is used. 
  
  Default: \code{NULL}.}
  %%%%%%%%%%%%%%%%%%%%%%%%%
  \item{timelab}{Time axis label. 
  
  Default: \code{"index"}; in case of a calendar axis: \code{"calendar date"}.}
  %%%%%%%%%%%%%%%%%%%%%%%%%
  \item{timetck}{length of tick marks on the time axis as a fraction of the smaller of the width or height of the plotting region; see \code{par}.  
  If \code{timetck >= 0.5}, \code{timetck} is interpreted as a fraction of the length of the time axis, so if \code{timetck = 1} 
  (and \code{timetcl = NULL}), vertical grid lines will be drawn. \cr
  Setting \code{timetck = NA} is to use \code{timetcl = -0.5} (which is the R default setting of \code{tck} and \code{tcl}).
  
  Default here: \code{0.02}.}
  \item{timetcl}{length of tick marks on the time axis as a fraction of the height of a line of text; see \code{par}. 
  With \code{timetcl = -0.5} (which is the R default setting of \code{tcl}), ticks will be drawn outward. 
  
  Default here: \code{0.5}.}
  %%%%%%%%%%%%%%%%%%%%%%%%%
  \item{spec.time.axis}{a list of tick mark and label specifications for individualized time axis labeling 
  (only effective if \code{label.time.axis = TRUE}):
  
  \itemize{
  \item[\code{at}:] locations of tick marks (when \code{NULL}, default plotting will be applied).
  Valid tick marks can be provided as numerical values or as dates. Dates are used only in the case \code{show.date = TRUE}, however, 
  and date formats should conform to \code{as.Date} or the format given in \code{date.format}. \cr
  Default: \code{NULL}. 
  \item[\code{labels}:] either a logical value specifying whether annotations at the tick marks are the tick marks themselves, 
  or any vector of labels. If \code{labels} is non-logical, \code{at} should be of same length. \cr
  Default: \code{TRUE}. 
  \item[\code{las}:] the style of axis labels, see \code{par}. \cr
  Default: \code{1} (always horizontal). 
  \item[\code{hadj}:] adjustment of labels horizontal to the reading direction, see \code{axis}. \cr
  Default: \code{NA} (centering is used). 
  \item[\code{padj}:] adjustment of labels perpendicular to the reading direction (this can be a vector of adjustments for each label), 
  see \code{axis}. \cr
  Default: \code{NA} (centering is used).                           
          }
  Mismatches will result in a reset to default plotting.                               
  } 
  %%%%%%%%%%%%%%%%%%%%%%%%%
  \item{label.period.axis}{Label the (Fourier) period axis? Logical. 
  
  Default: \code{TRUE}.}
  %%%%%%%%%%%%%%%%%%%%%%%%%
  \item{periodlab}{(Fourier) period axis label. 
  
  Default: \code{"period"}.} 
  %%%%%%%%%%%%%%%%%%%%%%%%%
  \item{periodtck}{length of tick marks on the period axis as a fraction of the smaller of the width or height 
  of the plotting region; see \code{par}. If \code{periodtck >= 0.5}, \code{periodtck} is interpreted as a fraction
  of the length of the period axis, so if \code{periodtck = 1} (and \code{periodtcl = NULL}), horizontal grid lines will be drawn.\cr 
  Setting \code{periodtck = NA} is to use \code{periodtcl = -0.5} (which is the R default setting of \code{tck} and \code{tcl}).
  
  Default here: \code{0.02}.}
  \item{periodtcl}{length of tick marks on the period axis as a fraction of the height of a line of text; see \code{par}. 
  With \code{periodtcl = -0.5} (which is the R default setting of \code{tcl}) ticks will be drawn outward. 
  
  Default here: \code{0.5}.}
  %%%%%%%%%%%%%%%%%%%%%%%%%  
  \item{spec.period.axis}{a list of tick mark and label specifications for individualized period axis labeling 
  (only effective if \code{label.period.axis = TRUE}):
  
  \itemize{
  \item[\code{at}:] locations of tick marks (when \code{NULL}, default plotting will be applied). Valid tick marks can be provided as 
  numerical and positive values only. \cr
  Default: \code{NULL}.        
  \item[\code{labels}:] either a logical value specifying whether annotations at the tick marks are the tick marks themselves, 
  or any vector of labels. If \code{labels} is non-logical, \code{at} should be of same length. \cr
  Default: \code{TRUE}. 
  \item[\code{las}:] the style of axis labels, see \code{par}. \cr
  Default: \code{1} (always horizontal). 
  \item[\code{hadj}:] adjustment of labels horizontal to the reading direction, see \code{axis}. \cr
  Default: \code{NA} (centering is used). 
  \item[\code{padj}:] adjustment of labels perpendicular to the reading direction (this can be a vector of adjustments for each label), 
  see \code{axis}. \cr
  Default: \code{NA} (centering is used).                           
          }  
  Mismatches will result in a reset to default plotting.                               
  }  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  \item{main}{an overall title for the plot. 
  
  Default: \code{NULL}.}
  \item{lwd}{line width of contour lines and ridge. 
  
  Default: \code{2}.}
  \item{lwd.axis}{line width of axes (image and legend bar). 
  
  Default: \code{1}.}
  \item{graphics.reset}{Reset graphical parameters? Logical. 
  
  Default: \code{TRUE}}
  \item{verbose}{Print verbose output on the screen? Logical. 
  
  Default: \code{FALSE}.}
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

Tian, H., and Cazelles, B., 2012. \code{WaveletCo}. 
Available at \url{https://cran.r-project.org/src/contrib/Archive/WaveletCo/}, archived April 2013; accessed July 26, 2013.

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
\code{\link{analyze.coherency}}, \code{\link{wc.avg}}, \code{\link{wc.sel.phases}}, \code{\link{wc.phasediff.image}}, 
\code{\link{wt.image}}, \code{\link{wt.avg}}, \cr
\code{\link{wt.sel.phases}}, \code{\link{wt.phase.image}}, \code{\link{reconstruct}}
}
\examples{
\dontrun{
## The following example is modified from Veleda et al., 2012:

series.length <- 3*128*24
x1 <- periodic.series(start.period = 1*24, length = series.length)
x2 <- periodic.series(start.period = 2*24, length = series.length)
x3 <- periodic.series(start.period = 4*24, length = series.length)
x4 <- periodic.series(start.period = 8*24, length = series.length)
x5 <- periodic.series(start.period = 16*24, length = series.length)
x6 <- periodic.series(start.period = 32*24, length = series.length)
x7 <- periodic.series(start.period = 64*24, length = series.length)
x8 <- periodic.series(start.period = 128*24, length = series.length)

x <- x1 + x2 + x3 + x4 + 3*x5 + x6 + x7 + x8 + rnorm(series.length)
y <- x1 + x2 + x3 + x4 - 3*x5 + x6 + 3*x7 + x8 + rnorm(series.length)

matplot(data.frame(x, y), type = "l", lty = 1, xaxs = "i", col = 1:2, 
 xlab = "index", ylab = "",
 main = "hourly series with periods of 1, 2, 4, 8, 16, 32, 64, 128 days", 
 sub = "(out of phase at period 16, different amplitudes at period 64)")
legend("topright", legend = c("x","y"), col = 1:2, lty = 1)

## The following dates refer to the local time zone 
## (possibly allowing for daylight saving time):      
my.date <- seq(as.POSIXct("2014-10-14 00:00:00", format = "\%F \%T"), 
               by = "hour", 
               length.out = series.length)     
my.data <- data.frame(date = my.date, x = x, y = y)

## Computation of cross-wavelet power and wavelet coherence, x over y:
## a natural choice of 'dt' in the case of hourly data is 'dt = 1/24',
## resulting in one time unit equaling one day. 
## This is also the time unit in which periods are measured.
my.wc <- analyze.coherency(my.data, c("x","y"), 
                           loess.span = 0, 
                           dt = 1/24, dj = 1/20, 
                           window.size.t = 1, window.size.s = 1/2, 
                           lowerPeriod = 1/4,
                           make.pval = TRUE, n.sim = 10)

## Plot of cross-wavelet power spectrum, 
## with color breakpoints according to quantiles:
wc.image(my.wc, 
   main = "cross-wavelet power spectrum, x over y",
   legend.params = list(lab = "cross-wavelet power levels (quantiles)"),
   periodlab = "period (days)")
## Note:
## The default time axis shows an index of given points in time, 
## which is the count of hours in our example.      
## By default, arrows are plotted which show the phase differences 
## of x over y at respective significant periods. 
## (Please see our guide booklet for further explanation.)

## The same plot, but with equidistant color breakpoints: 
wc.image(my.wc, color.key = "i", 
   main = "cross-wavelet power spectrum, x over y",
   legend.params = list(lab = "cross-wavelet power levels (equidistant)"),
   periodlab = "period (days)")

## The same plot, but adopting a palette of gray colors, 
## omitting the arrows:
wc.image(my.wc, color.key = "i", 
   main = "cross-wavelet power spectrum, x over y", 
   legend.params = list(lab = "cross-wavelet power levels (equidistant)"),
   color.palette = "gray( (1:n.levels)/n.levels )", 
   plot.arrow = FALSE,
   periodlab = "period (days)")
         
## The same plot, now with ridge of power:
wc.image(my.wc, color.key = "i", 
   main = "cross-wavelet power spectrum, x over y", 
   legend.params = list(lab = "cross-wavelet power levels (equidistant)"),
   color.palette = "gray( (1:n.levels)/n.levels )", 
   plot.arrow = FALSE,
   plot.ridge = TRUE, col.ridge = "red",
   periodlab = "period (days)")             
         
## The plot, turning back to arrows, now in yellow color:
wc.image(my.wc, color.key = "i", 
   main = "cross-wavelet power spectrum, x over y", 
   legend.params = list(lab = "cross-wavelet power levels (equidistant)"),
   color.palette = "gray( (1:n.levels)/n.levels )", 
   col.arrow = "yellow",
   periodlab = "period (days)")   
         
## Alternate styles of the time axis:          
    
## The plot with time elapsed in days, starting from 0 and proceeding 
## in steps of 50 days (50*24 hours), instead of the (default) time index:
index.ticks  <- seq(1, series.length, by = 50*24)
index.labels <- (index.ticks-1)/24    

## Insert your specification of time axis: 
wc.image(my.wc, color.key = "i", 
   main = "cross-wavelet power spectrum, x over y", 
   legend.params = list(lab = "cross-wavelet power levels (equidistant)"),
   color.palette = "gray( (1:n.levels)/n.levels )", 
   col.arrow = "yellow",
   periodlab = "period (days)", timelab = "time elapsed (days)",
   spec.time.axis = list(at = index.ticks, labels = index.labels))

## The plot with (automatically produced) calendar axis:
wc.image(my.wc, color.key = "i", 
   main = "cross-wavelet power spectrum, x over y", 
   legend.params = list(lab = "cross-wavelet power levels (equidistant)"),
   color.palette = "gray( (1:n.levels)/n.levels )", 
   col.arrow = "yellow",
   periodlab = "period (days)",
   show.date = TRUE, date.format = "\%F \%T")

## Individualizing your calendar axis (works with show.date = TRUE)...
## How to obtain, for example, monthly date ticks and labels:

## The sequence of tick positions:
monthly.ticks <- seq(as.POSIXct("2014-11-01 00:00:00", format = "\%F \%T"), 
                     as.POSIXct("2015-11-01 00:00:00", format = "\%F \%T"), 
                     by = "month")
## Observe that the following specification may produce an error:
## 'seq(as.Date("2014-11-01"), as.Date("2015-11-01"), by = "month")'
## Time of the day is missing here!

## The sequence of labels (e.g. information on month and year only):
monthly.labels <- strftime(monthly.ticks, format = "\%b \%Y")

## Insert your specification of time axis as parameter to wc.image: 
wc.image(my.wc, color.key = "i", 
   main = "cross-wavelet power spectrum, x over y", 
   legend.params = list(lab = "cross-wavelet power levels (equidistant)"),
   color.palette = "gray( (1:n.levels)/n.levels )", 
   col.arrow = "yellow", 
   periodlab = "period (days)",
   show.date = TRUE, date.format = "\%F \%T", 
   spec.time.axis = list(at = monthly.ticks, labels = monthly.labels, 
                         las = 2))
## Note: 
## The monthly ticks specify the midpoints of the colored cells and 
## match the location of corresponding (default) time index ticks.

## A cross-wavelet power plot with individualized period axis and exponent 
## to accentuate contrast in the image:
wc.image(my.wc, exponent = 0.5, color.key = "i", 
   main = "cross-wavelet power spectrum, x over y",
   legend.params = list(lab = "cross-wavelet power levels 
   (raised by exponent 0.5, equidistant levels)"),
   color.palette = "gray( (1:n.levels)/n.levels )", 
   col.arrow = "yellow", 
   periodlab = "period (days)",
   spec.period.axis = list(at = c(1,2,4,8,16,32,64,128))) 
         
## An option to switch to the corresponding frequency axis:
my.periods <- c(1,2,4,8,16,32,64,128)
my.frequencies <- paste("1/",my.periods, sep = "")
wc.image(my.wc, exponent = 0.5, color.key = "i", 
   main = "cross-wavelet power spectrum, x over y",
   legend.params = list(lab = "cross-wavelet power levels
   (raised by exponent 0.5, equidistant levels)"),
   color.palette = "gray( (1:n.levels)/n.levels )", 
   col.arrow = "yellow", 
   periodlab = "frequency (per day)", 
   spec.period.axis = list(at = my.periods, labels = my.frequencies))            

## Adding, for example, horizontal lines at period ticks...

## There is an option to add further objects to the image plot region, 
## by setting 'graphics.reset = FALSE' 
## (but recall previous par settings after plotting):

op <- par(no.readonly = TRUE)
wc.image(my.wc, exponent = 0.5, color.key = "i", 
   main = "cross-wavelet power spectrum, x over y",
   legend.params = list(lab = "cross-wavelet power levels
   (raised by exponent 0.5, equidistant levels)"),
   color.palette="gray( (1:n.levels)/n.levels )", 
   col.arrow = "yellow", 
   periodlab = "frequency (per day)", 
   spec.period.axis = list(at = my.periods, labels = my.frequencies),
   timelab = "",
   show.date = TRUE, date.format = "\%F \%T",
   graphics.reset = FALSE)            
abline(h = log2(my.periods))
year2015 <- as.POSIXct("2015-01-01 00:00:00", format = "\%F \%T")
abline(v = year2015)
axis(1, at = year2015, labels = 2015, padj = 1)
par(op)

## For further axis plotting options: 
## Please see the examples in our guide booklet,
## URL http://www.hs-stat.com/projects/WaveletComp/WaveletComp_guided_tour.pdf.    
            
## Plot of wavelet coherence of x over y, 
## with color breakpoints according to quantiles:
wc.image(my.wc, which.image = "wc", 
   main = "wavelet coherence, x over y",
   legend.params = list(lab = "wavelet coherence levels (quantiles)", 
                        lab.line = 3.5, label.digits = 3),
   periodlab = "period (days)")

## Plot of wavelet coherence, but with equidistant color breakpoints:
wc.image(my.wc, which.image = "wc", color.key = "i",
   main = "wavelet coherence, x over y",
   legend.params = list(lab = "wavelet coherence levels (equidistant)"),
   periodlab = "period (days)")  
   
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }