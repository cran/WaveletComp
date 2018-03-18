\name{reconstruct}
\alias{reconstruct}

\title{
Reconstruction of a (detrended) time series from output provided by an
object of class \code{"analyze.wavelet"} or \code{"analyze.coherency"}
}
\description{
This function reconstructs a (detrended) time series analyzed
by wavelet transformation using either function
\code{analyze.wavelet} or function \code{analyze.coherency}, subject to
optional criteria concerning: minimum wavelet power, significance of wavelet power at a
given significance level, specification of (Fourier) periods or
period bands, exclusive use of the power ridge and/or the cone of influence.
An option is provided to prevent the reconstructed series from final rescaling 
(applying the original (detrended) series' mean and standard deviation).

(If the object provided as input is of class \code{"analyze.coherency"}, 
then the number or name of the time series must be specified.)  

Optional: plot of wavelets used for reconstruction, plot of 
reconstructed series against original (detrended)
series. An option is given to individualize the time axis 
by specifying tick marks and labels.

Output includes the original (detrended) and the reconstructed
time series, along with reconstruction waves and parameters.
}

\usage{
reconstruct(WT, my.series = 1, lvl = 0, 
      only.coi = FALSE, 
      only.sig = TRUE, siglvl = 0.05, 
      only.ridge = FALSE, 
      sel.period = NULL, sel.lower = NULL, sel.upper = NULL,  
      rescale = TRUE,
      plot.waves = FALSE, plot.rec = TRUE, 
      lty = 1, lwd = 1, col = 1:2, ylim = NULL,
      show.legend = TRUE, 
      legend.coords = "topleft", legend.horiz = FALSE, legend.text = NULL,
      label.time.axis = TRUE, 
      show.date = FALSE, date.format = NULL, date.tz = NULL,
      timelab = NULL, timetck = 0.02, timetcl = 0.5,
      spec.time.axis = list(at = NULL, labels = TRUE, 
                            las = 1, hadj = NA, padj = NA),
      main.waves = NULL, main.rec = NULL, main = NULL, 
      lwd.axis = 1, 
      verbose = TRUE)
}

\arguments{
  \item{WT}{an object of class \code{"analyze.wavelet"} or \code{"analyze.coherency"}}
  \item{my.series}{In case \code{class(WT) = "analyze.coherency"}: number (\code{1} or \code{2}) or name of the series to be analyzed. 
  
  Default: \code{1}. }
  %%%%%%%%%%%%%%%%%%%%%%
  \item{lvl}{minimum level of wavelet power to be applied for the inclusion of reconstruction waves. 
  
  Default: \code{0}.}
  \item{only.coi}{Exclude borders influenced by edge effects in reconstruction, i.e. include the cone of influence only? Logical. 
  
  Default: \code{FALSE}.}
  \item{only.sig}{Use wavelet power significance in reconstruction? Logical. 
  
  Default: \code{TRUE}.}
  \item{siglvl}{level of wavelet power significance to be applied for the inclusion of reconstruction waves. 
  
  Default: 0.05.} 
  \item{only.ridge}{Select only the wavelet power ridge? Logical.
  
  Default: \code{FALSE}.}
  %%%%%%%%%%%%%%%%%%%%%%%%
  \item{sel.period}{a vector of numbers to select Fourier periods (or closest available periods)
  and corresponding wavelets for the reconstruction. 
  
  Default: \code{NULL}.}
  \item{sel.lower}{a number to define a lower Fourier period (or the closest available) for the selection 
  of a band of wavelets for the reconstruction.\cr 
  (Only effective if \code{sel.period = NULL}.) 
  
  Default: \code{NULL}.}
  \item{sel.upper}{a number to define an upper Fourier period (or the closest available) for the selection
  of a band of wavelets for the reconstruction.\cr 
  (Only effective if \code{sel.period = NULL}.) 
  
  Default: \code{NULL}.}
  %%%%%%%%%%%%%%%%%%%%%%%%
  \item{rescale}{Shall the reconstructed series finally be rescaled to attain the original (detrended) series' mean and standard deviation? Logical.
  
  Default: \code{TRUE}.}
  %%%%%%%%%%%%%%%%%%%%%%%%
  \item{plot.waves}{Shall reconstruction waves be plotted? Logical. 
  
  Default: \code{FALSE}.}
  \item{plot.rec}{Shall the reconstructed series (together with the original (detrended) series) be plotted? Logical. 
  
  Default: \code{TRUE}.}
  %%%%%%%%%%%%%%%%%%%%%%%%%
  \item{lty}{parameter for the plot of original vs. reconstructed series: line type, e.g. \code{1:2}. 
  
  Default: \code{1}.}
  \item{lwd}{parameter for the plot of original vs. reconstructed series: line width, e.g. \code{1:2}. 
  
  Default: \code{1}.}
  \item{col}{parameter for the plot of original vs. reconstructed series: color of lines. 
  
  Default: \code{1:2}.}
  %%%%%%%%%%%%%%%%%%%%%%%%%
  \item{ylim}{numeric vector of length \code{2}, providing the range of vertical coordinates for the plot. 
  
  Default: \code{NULL}.}
  %%%%%%%%%%%%%%%%%%%%%%%%%
  \item{show.legend}{Include legend into the plot of original vs. reconstructed series? Logical.
  
  Default: \code{TRUE}.}
  \item{legend.coords}{coordinates to position the legend (as in function \code{legend}). 
  
  Default: \code{"topleft"}.}
  \item{legend.horiz}{Set the legend horizontally rather than vertically? Logical. 
  
  Default: \code{FALSE}.}
  \item{legend.text}{legend text. Default: \code{c("original (detrended)", "reconstructed")}.}
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
  \item{main.waves}{an overall title for the plot of reconstruction waves. 
  
  Default: \code{NULL}.}
  \item{main.rec}{an overall title for the plot of original vs. reconstructed series. 
  
  Default: \code{NULL}.}
  \item{main}{an overall title for both plots. 
  
  Default: \code{NULL}.}
  \item{lwd.axis}{line width of axes. 
  
  Default: \code{1}.}
  \item{verbose}{Print verbose output on the screen? Logical. 
  
  Default: \code{TRUE}.}
}
% \details{ }
\value{A list of class \code{reconstruct} with the following elements:
   %%%%%%%%%%%%%%%%%
   
   \item{series}{a data frame building on \code{WT$series} with the following columns:
            \tabular{rlll}{
                     \tab date             \tab : \tab the calendar date (if available as column\cr
                     \tab                  \tab   \tab in WT$series) \cr
                     \tab <x>              \tab : \tab series <x>, with original name retained \cr
                     \tab                  \tab : \tab (detrended, if \code{loess.span != 0}) \cr
                     \tab <x>.trend        \tab : \tab the trend series (if \code{loess.span != 0}) \cr
                     \tab <x>.r            \tab : \tab the reconstructed (detrended) series
                     } 
                     
   Row names are taken over from WT$series, and so are dates if given as row names.
   If \code{WT} is of class \code{analyze.coherency}, the second series in the coherency analysis is retained; 
   if \code{loess.span != 0}, the second series is retained in the detrended version, and the trend is retained as well.   
   }
   %%%%%%%%%%%%%%%%%
   \item{rec.waves}{data frame of scaled waves used for reconstruction}
   %%%%%%%%%%%%%%%%%
   \item{loess.span}{parameter \code{alpha} in \code{loess} controlling the degree of time series smoothing, 
   if the time series was detrended; no detrending if \code{loess.span = 0}.}
   %%%%%%%%%%%%%%%%%
   \item{lvl}{minimum level of wavelet power for waves (wave segments) to be included
   in the reconstruction}
   \item{only.coi}{Was the influence of edge effects excluded? I.e. was the cone of influence used only?}  
   \item{only.sig}{Was wavelet power significance used in reconstruction?} 
   \item{siglvl}{level of wavelet power significance}
   \item{only.ridge}{Was the wavelet power ridge used only?}  
   %%%%%%%%%%%%%%%%%
   \item{rnum.used}{the vector of Fourier period numbers used for reconstruction }  
   %%%%%%%%%%%%%%%%%
   \item{rescale}{Was the reconstructed series rescaled according to the mean and standard deviation 
   taken from the original (detrended) series?} 
   %%%%%%%%%%%%%%%%%
   \item{dt}{time resolution, i.e. sampling resolution in the time domain, \code{1/dt} = number of observations per time unit}
   \item{dj}{frequency resolution, i.e. sampling resolution in the frequency domain, \code{1/dj} = number of suboctaves (voices per octave)}
   %%%%%%%%%%%%%%%%%
   \item{Period}{the Fourier periods 
   (measured in time units determined by \code{dt}, see the explanations concerning \code{dt})}
   \item{Scale}{the scales (the Fourier periods divided by the Fourier factor)}
   %%%%%%%%%%%%%%%%%
   \item{nc}{number of columns = number of observations = number of observation epochs; "epoch" meaning point in time}
   \item{nr}{number of rows = number of scales (Fourier periods)}
   %%%%%%%%%%%%%%%%%
   \item{axis.1}{tick levels corresponding to the time steps used for (cross-)wavelet transformation: \code{1, 1+dt, 1+2dt, ...}.
   The default time axis in plot functions provided by \code{WaveletComp} is determined by observation epochs, however; "epoch" meaning point in time. }
   \item{axis.2}{tick levels corresponding to the log of Fourier periods: \code{log2(Period)}. This determines the period axis in plot functions provided by \code{WaveletComp}.}
   %%%%%%%%%%%%%%%%%
   \item{date.format}{the format of calendar date (if available)}
   \item{date.tz}{the time zone of calendar date (if available)}
}
\references{
Carmona R., Hwang W.-L., and Torresani B., 1998. 
Practical Time Frequency Analysis. Gabor and Wavelet Transforms with an Implementation in S.
Academic Press, San Diego.

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
\code{\link{analyze.wavelet}}, \code{\link{wt.image}}, \code{\link{wt.avg}}, \code{\link{wt.sel.phases}}, \code{\link{wt.phase.image}}, \code{\link{analyze.coherency}}, 
\code{\link{wc.image}}, \code{\link{wc.avg}}, \code{\link{wc.sel.phases}}, \code{\link{wc.phasediff.image}}
}
\examples{
\dontrun{
## The following example is adopted from Liu et al., 2007:

series.length = 6*128*24
x1 <- periodic.series(start.period = 1*24, length = series.length)
x2 <- periodic.series(start.period = 8*24, length = series.length)
x3 <- periodic.series(start.period = 32*24, length = series.length)
x4 <- periodic.series(start.period = 128*24, length = series.length)

x <- x1 + x2 + x3 + x4

plot(x, type = "l", xlab = "index", ylab = "", xaxs = "i",
     main = "hourly series with periods of 1, 8, 32, 128 days")
  
my.data <- data.frame(x = x)

## Computation of wavelet power:
## a natural choice of 'dt' in the case of hourly data is 'dt = 1/24', 
## resulting in one time unit equaling one day. 
## This is also the time unit in which periods are measured.
my.w <- analyze.wavelet(my.data, "x", 
                       loess.span = 0, 
                       dt = 1/24, dj = 1/20, 
                       lowerPeriod = 1/4, 
                       make.pval = TRUE, n.sim = 10)

## Plot of wavelet power spectrum (with equidistant color breakpoints):  
wt.image(my.w, color.key = "interval", 
   legend.params = list(lab = "wavelet power levels"),
   periodlab = "period (days)")

## Reconstruction of the time series, 
## including significant components only: 
reconstruct(my.w)

## The same reconstruction, but showing wave components first:
reconstruct(my.w, plot.waves = TRUE)

## Reconstruction, including all components whether significant or not: 
reconstruct(my.w, only.sig = FALSE)

## Reconstruction, including significant components, 
## but selected periods only (e.g. ignoring period 8):  
reconstruct(my.w, sel.period = c(1,32,128))

## Reconstruction, including significant components,
## but the ridge only:
reconstruct(my.w, only.ridge = TRUE)

## Alternate styles of the time axis:

## The plot with time elapsed in days, starting from 0 and proceeding 
## in steps of 50 days (50*24 hours),
## instead of the (default) time index:
index.ticks  <- seq(1, series.length, by = 50*24)
index.labels <- (index.ticks-1)/24

## Insert your specification of time axis: 
reconstruct(my.w, only.ridge = TRUE,
   timelab = "time elapsed (days)",
   spec.time.axis = list(at = index.ticks, labels = index.labels))
            
## See the periods involved:
my.rec <- reconstruct(my.w, only.ridge = TRUE)
print(my.rec$Period[my.rec$rnum.used])

## The original and reconstructed time series can be retrieved:
plot(my.rec$series$x, type = "l", xlab = "index", ylab = "")
lines(my.rec$series$x.r, col="red")
legend("topleft", legend = c("original","reconstructed"), 
       lty = 1, col = c("black","red"))

## Please see also the examples in our guide booklet,
## URL http://www.hs-stat.com/projects/WaveletComp/WaveletComp_guided_tour.pdf.

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }