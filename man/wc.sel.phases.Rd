\name{wc.sel.phases}
\alias{wc.sel.phases}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Comparison plot of phases for selected periodic components of two time
series
}
\description{
This function plots the phases for periodic components of two time series, 
which are provided by an object of class \code{"analyze.coherency"}.

Periodic components can be selected by specification of a single
Fourier period or of a period band. In the latter case, and in the
default case (no specification at all), phases are averaged across
periods for each time series. Other options: restriction to the 
cone of influence, restriction to an area of significance (with respect 
to cross-wavelet power, wavelet coherence or individual wavelet power). 
Phase differences (i.e. angles, smoothed or not
smoothed) can be added to the plot.

(The time axis can be altered to display dates, see e.g. \code{wt.image}.
In particular, an option is given to individualize the phase and/or time axis
by specifying tick marks and labels.)
}
\usage{
wc.sel.phases(WC, sel.period = NULL, sel.lower = NULL, sel.upper = NULL, 
      only.coi = FALSE, 
      only.sig = TRUE, which.sig = "wp", siglvl = 0.05, 
      phase.cols = c("red", "blue"), 
      show.Angle = TRUE, use.sAngle = FALSE, Angle.col = "black", 
      show.legend = TRUE, legend.coords = "topleft", legend.horiz = TRUE,
      label.time.axis = TRUE, 
      show.date = FALSE, date.format = NULL, date.tz = NULL, 
      timelab = NULL, timetck = 0.02, timetcl = 0.5,
      spec.time.axis = list(at = NULL, labels = TRUE, 
                            las = 1, hadj = NA, padj = NA),
      label.phase.axis = TRUE, 
      phaselab = NULL, phasetck = 0.02, phasetcl = 0.5,
      spec.phase.axis = list(at = NULL, labels = TRUE, 
                             las = 1, hadj = NA, padj = NA),
      phaselim = c(-pi,pi+show.legend*ifelse(legend.horiz,0.8,2)),      
      main = NULL, sub = NULL,
      lwd = 1, lwd.Angle = 2, lwd.axis = 1, 
      verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{WC}{an object of class \code{"analyze.coherency"}. }
  %%%%%%%%%%%%%%%%%%%%%%
  \item{sel.period}{a single number which determines the (closest available) Fourier period to be selected. 
  
  Default: \code{NULL}.}
  \item{sel.lower}{a number to define a lower Fourier period (or the closest available) for the selection 
  of a band of periods (effective if \code{sel.period} is \code{NULL}). 
  
  Default: \code{NULL}.}
  \item{sel.upper}{a number to define an upper Fourier period (or the closest available) for the selection 
  of a band of periods (effective if \code{sel.period} is \code{NULL}). 
  
  Default: \code{NULL}.} 
  %%%%%%%%%%%%%%%%%%%%%%
  \item{only.coi}{Exclude borders influenced by edge effects, i.e. include the cone of influence only? Logical. 
  
  Default: \code{FALSE}.}
  \item{only.sig}{Use cross-wavelet power (or wavelet coherence, depending on \code{which.sig}) significance 
  to decide about the inclusion of (parts of) the series of phases? Logical. 
  
  Default: \code{TRUE}.}
  \item{which.sig}{Which spectrum should significance refer to?
                     
                       \tabular{rlll}{
                                \tab \code{"wp"} \tab : \tab cross-wavelet power \cr
                                \tab \code{"wc"} \tab : \tab wavelet coherence \cr
                                \tab \code{"wt"} \tab : \tab individual wavelet power
                                }
                                
  Default: \code{"wp"}
  }
  \item{siglvl}{level of cross-wavelet power (or wavelet coherence, depending on \code{which.sig}) significance.
  
  Default: \code{0.05}.} 
  %%%%%%%%%%%%%%%%%%%%%%
  \item{phase.cols}{a vector of two colors for the plot of (average) phases referring to the two time series.
  
  Default: \code{c("red","blue")}.}
  %%%%%%%%%%%%%%%%%%%%%%
  \item{show.Angle}{Show the (average) phase difference (the Angle) between the two series? Logical. 
  
  Default: \code{TRUE}.}
  \item{use.sAngle}{Use smoothed version of phase difference? Logical. 
  
  Default: \code{FALSE}.}
  \item{Angle.col}{color of the line of Angles. 
  
  Default: \code{"black"}.}
  %%%%%%%%%%%%%%%%%%%%%%%
  \item{show.legend}{Include legend? Logical.
  Default: \code{TRUE}.}
  \item{legend.coords}{Coordinates to position the legend (with the same options as given in function \code{legend}). 
  
  Default: \code{"topleft"}.}
  \item{legend.horiz}{Set the legend horizontally rather than vertically? Logical. 
  
  Default: \code{TRUE}.}
  %%%%%%%%%%%%%%%%%%%%%%%%
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
  \item{label.phase.axis}{Label the phase axis? Logical. 
  
  Default: \code{TRUE}.}
  \item{phaselab}{Phase axis label. 
  
  Default: \code{"phase"}.}
  %%%%%%%%%%%%%%%%%%%%%%%%%
  \item{phasetck}{length of tick marks on the phase axis as a fraction of the smaller of the width or height of the plotting region; see \code{par}.  
  If \code{phasetck >= 0.5}, \code{phasetck} is interpreted as a fraction of the length of the phase axis, so if \code{phasetck = 1} 
  (and \code{phasetcl = NULL}), horizontal grid lines will be drawn. \cr
  Setting \code{phasetck = NA} is to use \code{phasetcl = -0.5} (which is the R default setting of \code{tck} and \code{tcl}).
  
  Default here: \code{0.02}.}
  \item{phasetcl}{length of tick marks on the phase axis as a fraction of the height of a line of text; see \code{par}. 
  With \code{phasetcl = -0.5} (which is the R default setting of \code{tcl}), ticks will be drawn outward. 
  
  Default here: \code{0.5}.}
  %%%%%%%%%%%%%%%%%%%%%%%%%
  \item{spec.phase.axis}{a list of tick mark and label specifications for individualized phase axis labeling 
  (only effective if \code{label.phase.axis = TRUE}):
  
  \itemize{
  \item[\code{at}:] locations of tick marks (when \code{NULL}, default plotting will be applied). Valid tick marks can be provided as 
  numerical values between \code{-pi} and \code{pi}. \cr
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
  %%%%%%%%%%%%%%%%%%
  \item{phaselim}{numeric vector of length 2, giving the phase coordinate range.
  
  Default: \code{c(-pi,pi+0.8)}\cr
  (\code{+0.8} in order to accomodate the horizontal legend, \code{+2} in case of a vertical legend).}
  %%%%%%%%%%%%%%%%%%
  \item{main}{an overall title for the plot. 
  
  Default: \code{NULL}.}
  \item{sub}{a subtitle for the plot. 
  
  Default: \code{NULL}. In this case, the selected period range will be given in the subtitle.}
  \item{lwd}{width of lines of phases. 
  
  Default: \code{1}.}
  \item{lwd.Angle}{width of lines of (average) phase differences (the Angles) between the two series (this line will be plotted if \code{show.Angle = TRUE}).
  
  Default: \code{2}.}
  \item{lwd.axis}{line width of axes. 
  
  Default: \code{1}.}
  \item{verbose}{Print verbose output on the screen? Logical. 
  
  Default: \code{FALSE}.}
}
% \details{ }
\value{A list of class \code{"sel.phases"} with the following elements:
   \item{Period}{the selected period (or period band)}
   \item{Phase.x}{time series of (average) phases at the selected period (or period band), case of series x} 
   \item{Phase.y}{time series of (average) phases at the selected period (or period band), case of series y}
   \item{Angle}{time series of (average) phase differences (non-smoothed version) at the selected period (or period band)}
   \item{sAngle}{time series of (average) smoothed phase differences at the selected periods}
   \item{only.coi}{Is the influence of edge effects excluded? I.e. is the cone of influence used only?}  
   \item{only.sig}{Was significance used in selection of phases?} 
   \item{which.sig}{Which spectrum was used to refer to significance?  
                    
                       \tabular{rlll}{
                                \tab "wp" \tab : \tab cross-wavelet power \cr
                                \tab "wc" \tab : \tab wavelet coherence \cr
                                \tab "wt" \tab : \tab individual wavelet power
                                }
                                
   }
   \item{siglvl}{level of significance}
   \item{date}{time series of calendar date (if available)}
   \item{date.format}{the format of calendar date as provided}
   \item{date.tz}{the time zone of calendar date as provided}
   \item{axis.1}{tick levels corresponding to the time steps used for (cross-)wavelet transformation: 1, 1+dt, 1+2dt, ...}
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
\code{\link{analyze.coherency}}, \code{\link{wc.image}}, \code{\link{wc.avg}}, \code{\link{wc.phasediff.image}}, 
\code{\link{wt.image}}, \code{\link{wt.avg}}, \cr
\code{\link{wt.sel.phases}}, \code{\link{wt.phase.image}}, \code{\link{reconstruct}}
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

my.date <- seq(as.POSIXct("2014-10-14 00:00:00", format = "\%F \%T"), 
               by = "hour", 
               length.out = series.length)     
my.data <- data.frame(date = my.date, x = x, y = y)

## Computation of cross-wavelet power and wavelet coherency of x over y:
## a natural choice of 'dt' in the case of hourly data is 'dt = 1/24', 
## resulting in one time unit equaling one day.
## This is also the time unit in which periods are measured.
my.wc <- analyze.coherency(my.data, c("x","y"), loess.span = 0, 
                          dt = 1/24, dj = 1/20,
                          window.size.t = 1, window.size.s = 1/2, 
                          lowerPeriod = 1/4, 
                          make.pval = TRUE, n.sim = 10)

## Plot of cross-wavelet power spectrum, 
## with color breakpoints according to quantiles:
wc.image(my.wc, main = "cross-wavelet power spectrum, x over y",
   legend.params = list(lab = "cross-wavelet power levels (quantiles)"),
   periodlab = "period (days)")

## Select period 64 and compare plots of corresponding phases, including 
## the phase differences (angles) in their non-smoothed (default) version:
wc.sel.phases(my.wc, sel.period = 64, show.Angle = TRUE)

## With time elapsed in days
## (starting from 0 and proceeding in steps of 50 days) 
## instead of the (default) time index:
index.ticks  <- seq(1, series.length, by = 50*24)
index.labels <- (index.ticks-1)/24
wc.sel.phases(my.wc, sel.period = 64, show.Angle = TRUE,  
   timelab = "time elapsed (days)",
   spec.time.axis = list(at = index.ticks, labels = index.labels))

## The same plot, but with (automatically produced) calendar axis:
wc.sel.phases(my.wc, sel.period = 64, show.Angle = TRUE,   
   show.date = TRUE, date.format = "\%F \%T")

## For further axis plotting options:
## Please see the examples in our guide booklet,
## URL http://www.hs-stat.com/projects/WaveletComp/WaveletComp_guided_tour.pdf.    

## Now, select period 16...
## and observe that corresponding components are out of phase:
wc.sel.phases(my.wc, sel.period = 16, show.Angle = TRUE,
   show.date = TRUE, date.format = "\%F \%T")
## ... compare to period 4...
wc.sel.phases(my.wc, sel.period = 4, show.Angle = TRUE,
   show.date = TRUE, date.format = "\%F \%T")   

## In the following, no periods are selected. 
## In this case, instead of individual phases, the plot shows 
## average phases for each series: 
wc.sel.phases(my.wc)

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }