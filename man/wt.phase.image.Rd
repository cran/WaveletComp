\name{wt.phase.image}
\alias{wt.phase.image}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Image plot of the phases of periodic components for a single time
series
}
\description{
This function plots the wavelet phase image for a time series, which is provided by an object of
class \code{"analyze.wavelet"}, or alternatively of class \code{"analyze.coherency"}. 
(In the latter case, the series number or name must be specified.)  
The vertical axis shows the Fourier periods. The horizontal axis shows time step counts, but can
be easily transformed into a calendar axis if dates are provided in either row names or as a variable 
named \code{"date"} in the data frame at hand. Both axes can be relabeled.
In particular, an option is given to individualize the period and/or time axis 
by specifying tick marks and labels.

The color levels are defined according to equidistant breakpoints
(covering the interval from \code{-pi} to \code{+pi}), with the number of levels
as a further parameter. In addition, there is an option to adopt
an individual color palette.

If the default palette is retained, colors indicate the following. Green: Phases
close to zero. Red: phases close to \code{+pi}. Blue: phases close to \code{-pi}.

Further plot design options concern: plot of the cone of
influence, plot of contour lines to border areas of significance
with respect to cross-wavelet power or wavelet coherency at a
given significance level, plot of power ridges.

Finally, there is an option to insert and format a color legend (a
right-hand vertical color bar) and to set the plot title.  For
further processing of the plot, graphical parameters of plot
regions are provided as output.
}
\usage{
wt.phase.image(WT, my.series = 1,
       plot.coi = TRUE, 
       plot.contour = TRUE, 
       siglvl = 0.1, col.contour = "white",
       plot.ridge = TRUE, col.ridge = "black", 
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
  \item{WT}{an object of class \code{"analyze.wavelet"} or \code{"analyze.coherency"}}
  \item{my.series}{In case \code{class(WT) = "analyze.coherency"}: number (\code{1} or \code{2}) or name of the series to be analyzed. 
  
  Default: \code{1}. }
  %%%%%%%%%%%%%%%%%%%%%%
  \item{plot.coi}{Plot cone of influence? Logical. Default: \code{TRUE}}
  \item{plot.contour}{Plot contour lines to border the area of wavelet power significance? Logical.
  
  Default: \code{TRUE}.
  }
  \item{siglvl}{level of wavelet power significance to be applied to the plot of contour lines. 
  
  Default: \code{0.1}.
  }
  \item{col.contour}{color of contour lines. Default: \code{"white"}.}
  %%%%%%%%%%%%%%%%%%%%%%%%%
  \item{plot.ridge}{Plot the wavelet power ridge? Logical.
  
  Default: \code{TRUE}.
  }
  \item{col.ridge}{ridge color. Default: \code{"black"}.}
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
  stored in \code{WT}. If unspecified, date formatting is attempted according to \code{as.Date}. 
  
  Default: \code{NULL}.}
  \item{date.tz}{a character string specifying the time zone of calendar date; see \code{strptime}. Explicit information given here will overturn 
  any specification stored in \code{WT}. If unspecified, \code{""} (the local time zone) is used. 
  
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

Torrence C., and Compo G.P., 1998. 
A practical guide to wavelet analysis. 
Bulletin of the American Meteorological Society 79 (1), 61--78.
}
\author{Angi Roesch and Harald Schmidbauer}
% \note{ }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ 
\code{\link{analyze.wavelet}}, \code{\link{wt.image}}, \code{\link{wt.avg}}, \code{\link{wt.sel.phases}}, \code{\link{reconstruct}} 
}
\examples{ 
\dontrun{
## The following example is adopted from Liu et al., 2007:

series.length <- 6*128*24
x1 <- periodic.series(start.period = 1*24, length = series.length)
x2 <- periodic.series(start.period = 8*24, length = series.length)
x3 <- periodic.series(start.period = 32*24, length = series.length)
x4 <- periodic.series(start.period = 128*24, length = series.length)

x <- x1 + x2 + x3 + x4

plot(x, type = "l", xlab = "index", ylab = "", xaxs = "i",
     main = "hourly series with periods of 1, 8, 32, 128 days")

my.date <- seq(as.POSIXct("2014-10-14 00:00:00", format = "\%F \%T"), 
               by = "hour", 
               length.out = series.length)     
my.data <- data.frame(date = my.date, x = x)  

## Computation of wavelet power:
## a natural choice of 'dt' in the case of hourly data is 'dt = 1/24', 
## resulting in one time unit equaling one day.
## This is also the time unit in which periods are measured.
my.wt <- analyze.wavelet(my.data, "x", 
                        loess.span = 0, 
                        dt = 1/24, dj = 1/20, 
                        lowerPeriod = 1/4, 
                        make.pval = TRUE, n.sim = 10)

## Plot of wavelet power spectrum with equidistant color breakpoints:  
wt.image(my.wt, color.key = "i", main = "wavelet power spectrum",
   legend.params = list(lab = "wavelet power levels (equidistant)"),
   periodlab = "period (days)")

## Default image of phases:
wt.phase.image(my.wt, 
   main = "image of phases", 
   periodlab = "period (days)")

## With time elapsed in days 
## (starting from 0 and proceeding in steps of 50 days) 
## instead of the (default) time index:
index.ticks  <- seq(1, series.length, by = 50*24)
index.labels <- (index.ticks-1)/24
wt.phase.image(my.wt, 
   main = "image of phases", 
   periodlab = "period (days)",
   timelab = "time elapsed (days)",
   spec.time.axis = list(at = index.ticks, labels = index.labels))

## The same plot, but with (automatically produced) calendar axis:
wt.phase.image(my.wt, 
   main = "image of phases", periodlab = "period (days)", 
   show.date = TRUE, date.format = "\%F \%T")

## For further axis plotting options:
## Please see the examples in our guide booklet,
## URL http://www.hs-stat.com/projects/WaveletComp/WaveletComp_guided_tour.pdf.    

## Image plot of phases with numerals as labels of the color legend bar: 
wt.phase.image(my.wt, 
   legend.params=list(pi.style = FALSE, label.digits = 2))
   
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }