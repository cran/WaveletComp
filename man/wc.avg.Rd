\name{wc.avg}
\alias{wc.avg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot cross-wavelet power averages and wavelet coherence averages across time of two time series
}
\description{
This function plots cross-wavelet power averages across time, or alternatively 
wavelet coherence averages, of two time series, 
which are provided by an object of class \code{"analyze.coherency"}. 
The vertical axis shows the Fourier periods. The horizontal axis shows the averages.
User-defined minimum and maximum levels can be applied to cross-wavelet power averages, 
minimum levels can be applied to coherence averages. 
Also, an option is given to individualize the period axis and/or axis of averages 
by specifying tick marks and labels.

There is an option to label periods according to significance of
averages (if p-values are provided) at given levels of significance. 
Labels are point symbols along the line of averages which can be assigned individually.

The idea to show significance levels by colors of plotting characters and its implementation 
has been adopted from Huidong Tian and Bernard Cazelles (archived R package \code{WaveletCo}).
}
\usage{
wc.avg(WC, which.avg = "wp", exponent = 1,
   show.siglvl = TRUE, 
   siglvl = c(0.05, 0.1), 
   sigcol = c("red", "blue"), sigpch = 20, sigcex = 1,
   minimum.level = NULL, maximum.level = NULL,
   label.avg.axis = TRUE, 
   averagelab = NULL, averagetck = 0.02, averagetcl = 0.5,
   spec.avg.axis = list(at = NULL, labels = TRUE, 
                        las = 1, hadj = NA, padj = NA),
   label.period.axis = TRUE, 
   periodlab = NULL, periodtck = 0.02, periodtcl = 0.5,
   spec.period.axis = list(at = NULL, labels = TRUE, 
                           las = 1, hadj = NA, padj = NA),
   show.legend = TRUE, legend.coords = "topright", 
   main = NULL, 
   lwd = 1, col = 1, 
   lwd.axis = 1,
   verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{WC}{an object of class \code{"analyze.coherency"}.}
  \item{which.avg}{Which averages should be plotted?
  
                   \tabular{rlll}{
                            \tab \code{"wp"} \tab : \tab cross-wavelet power \cr
                            \tab \code{"wc"} \tab : \tab wavelet coherence                            
                           }
                                
  Default: \code{"wp"}.                
  }
  %%%%%%%%%%%%%%%%%%%%%%%%%
  \item{exponent}{Exponent applied to averages before plotting; the exponent should be positive.
  
  Default: \code{1}.
  }
  %%%%%%%%%%%%%%%%%%%%%%%%%
  \item{show.siglvl}{Label periods according to significance of averages? (Effective only if p-values are provided.) 
  
  Default: \code{TRUE}.
  }
  \item{siglvl}{a vector of significance levels (of any length and order). Default: \code{c(0.05, 0.1)}.}
  \item{sigcol}{a vector of colors (should be of same length as and correspond to \code{siglvl}, 
  otherwise colors \code{1 : length(siglvl)}).
  
  Default: \code{c("red","blue")}.
  }
  \item{sigpch}{a vector of plotting "characters" (symbols) to use as labels of significance. 
  (It should be of same length as and correspond to \code{siglvl} to produce different plotting labels, 
  otherwise the default setting is used. A single input value affects all labels.) 
  
  Default: \code{20}.}
  \item{sigcex}{a numerical vector working as size of labels of significance. 
  (It should be of same length as and correspond to \code{siglvl} to produce different-sized labels, 
  otherwise the default setting is used.  A single input value affects all labels. 
  Note that \code{sigcex} is affected by \code{cex} in \code{par}.) 
  
  Default: \code{1}.}
  %%%%%%%%%%%%%%%%%%%%%%
  \item{minimum.level}{Minimum plot level of cross-wavelet power or wavelet coherence averages considered. 
  
  Default: NULL (referring to minimum level observed).}
  \item{maximum.level}{Maximum plot level of cross-wavelet power averages considered. 
  
  Default: NULL (referring to maximum level observed). Wavelet coherence has maximum average level \code{1} by definition.}
  %%%%%%%%%%%%%%%%%%%%%%
  \item{label.avg.axis}{Label the axis of averages? Logical. 
  
  Default: \code{TRUE}.}
  \item{averagelab}{Label for the axis of averages. 
  
  Default: \code{"average cross-wavelet power"} or \code{"average coherence"}, depending on \code{which.avg}.}
  %%%%%%%%%%%%%%%%%%%%%%
  \item{averagetck}{length of tick marks on the axis of averages as a fraction of the smaller of the width or height of 
  the plotting region; see \code{par}.  
  If \code{averagetck >= 0.5}, \code{averagetck} is interpreted as a fraction of the length of the axis of averages, 
  so if \code{averagetck = 1} (and \code{averagetcl = NULL}), vertical grid lines will be drawn. \cr
  Setting \code{averagetck = NA} is to use \code{averagetcl = -0.5} (which is the R default setting of \code{tck} and \code{tcl}).
  
  Default here: \code{0.02}.}
  \item{averagetcl}{length of tick marks on the axis of averages as a fraction of the height of a line of text; see \code{par}. 
  With \code{averagetcl = -0.5} (which is the R default setting of \code{tcl}) ticks will be drawn outward.
  
  Default here: \code{0.5}.}
  %%%%%%%%%%%%%%%%%%%%%%
  \item{spec.avg.axis}{a list of tick mark and label specifications for individualized labeling of the axis of averages 
  (only effective if \code{label.avg.axis = TRUE}):
  
  \itemize{
   \item[\code{at}:] locations of tick marks (when \code{NULL}, default plotting will be applied). Valid tick marks can be provided as 
   numerical and non-negative values only. \cr
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
  %%%%%%%%%%%%%%%%%%%%%%%
  \item{label.period.axis}{Label the (Fourier) period axis? Logical. 
  
  Default: \code{TRUE}.}
  \item{periodlab}{(Fourier) period axis label. 
  
  Default: \code{"period"}.} 
  %%%%%%%%%%%%%%%%%%%%%%%%
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
  %%%%%%%%%%%%%%%%%%%%%%%%%                             
  \item{show.legend}{Include legend of significance levels into the plot? Logical.
  
  Default: \code{TRUE}.}
  \item{legend.coords}{coordinates to position the legend (as in function \code{legend}). 
  
  Default: \code{"topright"}.}
  \item{main}{an overall title for the plot. 
  
  Default: NULL.}
  \item{lwd}{width of line of averages. 
  
  Default: 1.}
  \item{col}{color of line of averages. 
  
  Default: \code{"black"}.}
  \item{lwd.axis}{line width of axes. 
  
  Default: 1.}
  \item{verbose}{Print verbose output on the screen? Logical. 
  
  Default: \code{FALSE}.}
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

Tian, H., and Cazelles, B., 2012. \code{WaveletCo}. Available at \url{https://cran.r-project.org/src/contrib/Archive/WaveletCo/}, archived April 2013; accessed July 26, 2013.

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
\code{\link{analyze.coherency}}, \code{\link{wc.image}}, \code{\link{wc.sel.phases}}, \code{\link{wc.phasediff.image}}, \code{\link{wt.image}},
\code{\link{wt.avg}}, \cr \code{\link{wt.sel.phases}}, \code{\link{wt.phase.image}}, \code{\link{reconstruct}}
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

## Plot of cross-wavelet power, 
## with color breakpoints according to quantiles:
wc.image(my.wc, main = "cross-wavelet power spectrum, x over y",
   legend.params = list(lab = "cross-wavelet power levels (quantiles)"),
   periodlab = "period (days)")
## Note:
## The default time axis shows an index of given points in time, 
## which is the count of hours in our example.      
## By default, arrows are plotted which show the phase differences 
## of x over y at respective significant periods. 
## (Please see our guide booklet for further explanation.)
         
## With time elapsed in days 
## (starting from 0 and proceeding in steps of 50 days)
## instead of the (default) time index:
index.ticks  <- seq(1, series.length, by = 50*24)
index.labels <- (index.ticks-1)/24

## Insert your specification of the time axis: 
wc.image(my.wc, color.key = "i", 
   main = "cross-wavelet power spectrum, x over y", 
   legend.params = list(lab = "cross-wavelet power levels (quantiles)"),
   periodlab = "period (days)", timelab = "time elapsed (days)",
   spec.time.axis = list(at = index.ticks, labels = index.labels))
       
## Plot of average cross-wavelet power:
wc.avg(my.wc, siglvl = 0.05, sigcol = "red", periodlab = "period (days)")

## The same plot, but with enhanced symbol size, user-defined axes, 
## minimum and a maximum plot level of averages:
wc.avg(my.wc, siglvl = 0.05, sigcol = "red", sigcex = 1.5, 
   minimum.level = 0, maximum.level = 17, 
   periodlab = "period (days)",
   spec.period.axis = list(at = c(1,2,4,8,16,32,64,128)), 
   spec.avg.axis = list(at = seq(0,16,2)), 
   lwd = 1.5)
       
## Another style of the plot:
## 'cex.axis' in 'par' controls for the size of axis tick labels, 
## while 'cex.lab' controls for the size of axis and legend labels. 
## Note that scaling by 'cex' would also affect 'sigcex'. 
op <- par(no.readonly = TRUE)
par(cex.lab = 1.3, cex.axis = 1.1)
wc.avg(my.wc, siglvl = 0.05, sigcol = "red", sigcex = 1.5, 
   minimum.level = 0, maximum.level = 17, 
   periodlab = "period (days)",
   spec.period.axis = list(at = c(1,2,4,8,16,32,64,128)), 
   spec.avg.axis = list(at = seq(0,16,2)),
   lwd = 1.5)
par(op)      

## Plot of wavelet coherence 
## (with color breakpoints according to quantiles):
wc.image(my.wc, which.image = "wc", main = "wavelet coherence, x over y",
   legend.params = list(label.digits = 3),
   periodlab = "period (days)")
         
## Plot of average wavelet coherence:
wc.avg(my.wc, which.avg = "wc", 
   siglvl = 0.05, sigcol = "red", legend.coords = "topleft", 
   periodlab = "period (days)", 
   lwd = 1.5)

## The same plot, setting the minimum plot level of averages to 0: 
wc.avg(my.wc, which.avg = "wc", 
   siglvl = 0.05, sigcol = "red", legend.coords = "topleft", 
   minimum.level = 0, 
   periodlab = "period (days)", 
   lwd = 1.5)
 
## Please see also the examples in our guide booklet,
## URL http://www.hs-stat.com/projects/WaveletComp/WaveletComp_guided_tour.pdf.

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }