\name{wc.phasediff.image}
\alias{wc.phasediff.image}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Image plot of phase differences of periodic components for two time
series
}
\description{
This function plots the phase difference image of two time series, which is 
provided by an object of class \code{"analyze.coherency"}.  
The vertical axis shows the Fourier periods. The horizontal axis shows time step counts, but can 
be easily transformed into a calendar axis if dates are provided in either row names or a variable 
named \code{"date"} in the data frame at hand. Both axes can be relabeled.
In particular, an option is given to individualize the period and/or time axis 
by specifying tick marks and labels.

The color levels are defined according to equidistant breakpoints
(covering the interval from \code{-pi} to \code{+pi}), with the number of levels
as a further parameter. In addition, there is an option to adopt
an individual color palette.

If the default palette is retained, colors indicate the following. Green:
phase differences close to zero, which means that the two time
series are in phase at the respective period. Yellowgreen: in
phase, series 1 leading. Turquoise: in phase, series 2 leading.
Red: phase differences are close to \code{+pi}, out of phase, series
2 leading. Blue: phase differences are close to \code{-pi}, out of
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
wc.phasediff.image(WC, use.sAngle = FALSE, 
   plot.coi = TRUE, 
   plot.contour = TRUE, which.contour = "wp", 
   siglvl = 0.1, col.contour = "white",
   n.levels = 100, 
   color.palette = "rainbow(n.levels, start = 0, end = .7)", 
   useRaster = TRUE, max.contour.segments = 250000,
   plot.legend = TRUE,
   legend.params = list(width = 1.2, shrink = 0.9, mar = 5.1, 
                        n.ticks = 6, 
                        pi.style = TRUE, 
                        label.digits = 1, label.format = "f", 
                        lab = NULL, lab.line = 3),
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
  %%%%%%%%%%%%%%%%%%%
  \item{use.sAngle}{Use smoothed version of phase difference? Logical. 
  
  Default: \code{FALSE}.}
  %%%%%%%%%%%%%%%%%%%%%%%%%
  \item{plot.coi}{Plot cone of influence? Logical. Default: \code{TRUE}}
  %%%%%%%%%%%%%%%%%%%
  \item{plot.contour}{Plot contour lines to border the area of cross-wavelet power (or wavelet coherence, 
  depending on \code{which.contour}) significance? Logical.
  
  Default: \code{TRUE}.
  }
  \item{which.contour}{Contour lines of which spectrum should be plotted?
  
                       \tabular{rlll}{
                                \tab \code{"wp"} \tab : \tab cross-wavelet power \cr
                                \tab \code{"wc"} \tab : \tab wavelet coherence 
                                }
                                
  Default: \code{"wp"}.
                       }
  \item{siglvl}{level of cross-wavelet power (or wavelet coherence, depending on \code{which.contour}) 
  significance to be applied to the plot of contour lines. 
  
  Default: \code{0.1}.
  }
  \item{col.contour}{color of contour lines. Default: \code{"white"}.}
  %%%%%%%%%%%%%%%%%%%%%%%%%
   \item{n.levels}{Number of color levels. Default: \code{100}.}
  \item{color.palette}{Definition of color levels. (The color palette will be assigned to levels in reverse order!)
  
  Default: \code{"rainbow(n.levels, start = 0, end = .7)"}.
  }
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
  \item[\code{pi.style}:] style of labels, logical: if \code{TRUE}, the symbol \code{"pi"} is used to label the legend bar; otherwise, labels are numerals. \cr
  Default: \code{TRUE}. 
  \item[\code{label.digits}:] number of digits of (numerical factors of) labels. \cr
  Default: \code{1} if \code{pi.style} is \code{TRUE}, else \code{2}. 
  \item[\code{label.format}:] format of labels. \cr
  Default: \code{"f"}. 
  \item[\code{lab}:] axis label. \cr
  Default: \code{NULL}. 
  \item[\code{lab.line}:] line (in user coordinate units) where to put the axis label. \cr
  Default: \code{3}.
          }         
  }
  %%%%%%%%%%%%%%%%%%%%%%%%%
  \item{label.time.axis}{Label the time axis? Logical. 
  
  Default: \code{TRUE}.}
  \item{show.date}{Show calendar dates? (Effective only if dates are available as row names or by variable 
  \code{date} in the data frame which is analyzed.) Logical. 
  
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
  of the length of the period axis, so if \code{periodtck = 1} (and \code{periodtcl = NULL}), horizontal grid lines will be drawn. \cr
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
  %%%%%%%%%%%%%%%%%%%%%%%%%% 
  \item{main}{an overall title for the plot. 
  
  Default: \code{NULL}.}
  \item{lwd}{line width of contour lines and ridge. 
  
  Default: \code{2}.}
  \item{lwd.axis}{line width of axes (image and legend bar). 
  
  Default: \code{1}.}
  \item{graphics.reset}{Reset graphical parameters? Logical. 
  
  Default: \code{TRUE}.}
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
\code{\link{analyze.coherency}}, \code{\link{wc.image}}, \code{\link{wc.avg}}, \code{\link{wc.sel.phases}}, \code{\link{wt.image}}, \code{\link{wt.avg}}, \code{\link{wt.sel.phases}}, \code{\link{wt.phase.image}}, \code{\link{reconstruct}}
}
\examples{
\dontrun{
## The following example is modified from Veleda et al., 2012:

series.length <- 3*128*24
x1 <- periodic.series(start.period = 1*24, length = series.length)
x2 <- periodic.series(start.period = 2*24, length = series.length)
x3a <- periodic.series(start.period = 4*24, length = series.length)
x3b <- periodic.series(start.period = 4*24, length = series.length,
                       phase = 24)
x4 <- periodic.series(start.period = 8*24, length = series.length)
x5 <- periodic.series(start.period = 16*24, length = series.length)
x6 <- periodic.series(start.period = 32*24, length = series.length)
x7 <- periodic.series(start.period = 64*24, length = series.length)
x8 <- periodic.series(start.period = 128*24, length = series.length)

x <- x1 + x2 + x3a + x4 + 3*x5 + x6 + x7 + x8 + rnorm(series.length)
y <- x1 + x2 + x3b + x4 - 3*x5 + x6 + 3*x7 + x8 + rnorm(series.length)

matplot(data.frame(x, y), type = "l", lty = 1, xaxs = "i", col = 1:2, 
 xlab = "index", ylab = "",
 main = "hourly series with periods of 1, 2, 4, 8, 16, 32, 64, 128 days", 
 sub = "(different phases at periods 4 and 16)")
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

## Default plot of phase differences 
## (with contour lines referring to cross-wavelet power)
wc.phasediff.image(my.wc, which.contour = "wp", 
   main = "image of phase differences, x over y",
   periodlab = "period (days)")
      
## With time elapsed in days 
## (starting from 0 and proceeding in steps of 50 days) 
## instead of the (default) time index:
index.ticks  <- seq(1, series.length, by = 50*24)
index.labels <- (index.ticks-1)/24      
wc.phasediff.image(my.wc, which.contour = "wp", 
   main = "image of phase differences, x over y",
   periodlab = "period (days)",
   timelab = "time elapsed (days)",
   spec.time.axis = list(at = index.ticks, labels = index.labels))

## The same plot, but with (automatically produced) calendar axis:  
wc.phasediff.image(my.wc, which.contour = "wp", 
   main = "image of phase differences, x over y",
   periodlab = "period (days)",
   show.date = TRUE, date.format = "\%F \%T")
     
## For further axis plotting options:
## Please see the examples in our guide booklet,
## URL http://www.hs-stat.com/projects/WaveletComp/WaveletComp_guided_tour.pdf. 

## Plot of phase difference with numerals as labels of the color legend bar: 
wc.phasediff.image(my.wc, 
   legend.params=list(pi.style = FALSE, label.digits = 2))
   
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }