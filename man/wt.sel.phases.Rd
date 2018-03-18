\name{wt.sel.phases}
\alias{wt.sel.phases}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot phases for selected periodic components of a single time series
}
\description{
This function plots the phases for selected periodic components of a time series, 
which are provided by an object of class \code{"analyze.wavelet"}, or alternatively of class
\code{"analyze.coherency"}. (In the latter case, the series number or name must be specified.)

Periodic components can be selected by specification of a single
Fourier period or of a period band. In the latter case, and in the
default case (no specification at all), there is an option to
average the phases across periods. Other options: restriction to the 
cone of influence, restriction to an area of significance (with respect 
to wavelet power).

(The time axis can be altered to give dates, see e.g. \code{wt.image}.
In particular, an option is given to individualize the phase and/or time axis
by specifying tick marks and labels.)
}
\usage{
wt.sel.phases(WT, my.series = 1,
      sel.period = NULL, sel.lower = NULL, sel.upper = NULL, 
      only.coi = FALSE,
      only.sig = TRUE, siglvl = 0.05, 
      show.avg.phase = FALSE, phase.avg.col = "black", 
      label.time.axis = TRUE, 
      show.date = FALSE, date.format = NULL, date.tz = NULL, 
      timelab = NULL, timetck = 0.02, timetcl = 0.5,
      spec.time.axis = list(at = NULL, labels = TRUE, 
                            las = 1, hadj = NA, padj = NA),
      label.phase.axis = TRUE, 
      phaselab = NULL, phasetck = 0.02, phasetcl = 0.5,
      spec.phase.axis = list(at = NULL, labels = TRUE, 
                             las = 1, hadj = NA, padj = NA),
      main = NULL, sub = NULL, 
      lwd = 1, lwd.axis = 1,
      verbose = FALSE)
}
\arguments{
  \item{WT}{an object of class \code{"analyze.wavelet"} or \code{"analyze.coherency"}.}
  \item{my.series}{In case \code{class(WT) = "analyze.coherency"}: number (\code{1} or \code{2}) or name of the series to be analyzed. 
  
  Default: \code{1}. }
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
  \item{only.sig}{Use wavelet power significance to decide about the inclusion of (parts of) the phases' series? Logical. 
  
  Default: \code{TRUE}.}
  \item{siglvl}{level of wavelet power significance. Default: \code{0.05}.} 
  %%%%%%%%%%%%%%%%%%%%%%
  \item{show.avg.phase}{Show average phases over selected periods? (Effective only if a band of periods is selected.) Logical. 
  
  Default: \code{FALSE}.}
  \item{phase.avg.col}{color of line of phase averages. Default: "black".}
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
  \item{main}{an overall title for the plot. 
  
  Default: \code{NULL}.}
  \item{sub}{a subtitle for the plot. 
  
  Default: \code{NULL}. In this case, the selected period range will be given in the subtitle.}
  \item{lwd}{line width of phases. 
  
  Default: \code{1}.}
  \item{lwd.axis}{line width of axes. 
  
  Default: \code{1}.}
  \item{verbose}{Print verbose output on the screen? Logical.
  
  Default: \code{FALSE}.}
}
% \details{ }
\value{A list of class \code{"sel.phases"} with the following elements:
   %%%%%%%%%%%%%%%%%
   \item{Period}{the selected period (or period band)}
   \item{Phase}{time series of (average) phases at the selected period (or period band)} 
   \item{only.coi}{Is the influence of edge effects excluded? I.e. is the cone of influence used only?}  
   \item{only.sig}{Was wavelet power significance used in selection of phases?} 
   \item{siglvl}{level of wavelet power significance}
   %%%%%%%%%%%%%%%%%
   \item{date}{time series of calendar date (if available)}
   \item{date.format}{the format of calendar date (if available)}
   \item{date.tz}{the time zone of calendar date (if available)}
   %%%%%%%%%%%%%%%%%
   \item{axis.1}{tick levels corresponding to the time steps used for (cross-)wavelet transformation: \code{1, 1+dt, 1+2dt, ...}.
   The default time axis in plot functions provided by \code{WaveletComp} is determined by observation epochs, however; "epoch" meaning point in time. }
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
\code{\link{analyze.wavelet}}, \code{\link{wt.image}}, \code{\link{wt.avg}}, \code{\link{wt.phase.image}}, 
\code{\link{reconstruct}} 
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

## Select period 8 and plot corresponding phases across time:
wt.sel.phases(my.wt, sel.period = 8)

## With time elapsed in days 
## (starting from 0 and proceeding in steps of 50 days) 
## instead of the (default) time index:
index.ticks  <- seq(1, series.length, by = 50*24)
index.labels <- (index.ticks-1)/24
wt.sel.phases(my.wt, sel.period = 8,  
   timelab = "time elapsed (days)",
   spec.time.axis = list(at = index.ticks, labels = index.labels))

## The same plot, but with (automatically produced) calendar axis:
wt.sel.phases(my.wt, sel.period = 8,  
   show.date = TRUE, date.format = "\%F \%T")

## For further axis plotting options:
## Please see the examples in our guide booklet,
## URL http://www.hs-stat.com/projects/WaveletComp/WaveletComp_guided_tour.pdf.        

## In the following, no period is selected.
## By setting 'show.avg.phase = TRUE', the plot shows average phases 
## instead of individual phases: 
wt.sel.phases(my.wt, show.avg.phase = TRUE)

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}