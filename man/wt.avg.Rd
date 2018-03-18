\name{wt.avg}
\alias{wt.avg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of wavelet power averages across time of a single time series
}
\description{
This function plots wavelet power averages across time of a single time
series, which are provided by an object of class \code{"analyze.wavelet"}, 
or alternatively of class \code{"analyze.coherency"}. 
(In the latter case, the series number or name must be specified.) 
The vertical axis shows the Fourier periods. The horizontal axis shows the averages.
User-defined minimum and maximum average levels can be applied. 
Also, an option is given to individualize the period axis and/or axis of averages 
by specifying tick marks and labels.

There is an option to label periods according to significance of
averages (if p-values are provided) at given levels of significance.
Labels are point symbols along the
line of averages which can be assigned individually.

The idea to show significance levels by colors of plotting characters and its 
implementation has been adopted from Huidong Tian and Bernard Cazelles 
(archived R package \code{WaveletCo}).
}
\usage{
wt.avg(WT, my.series = 1, exponent = 1,
  show.siglvl = TRUE, siglvl = c(0.05, 0.1), 
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
  \item{WT}{an object of class \code{"analyze.wavelet"} or \code{"analyze.coherency"}}
  \item{my.series}{In case \code{class(WT) = "analyze.coherency"}: number (\code{1} or \code{2}) or name of 
  the series to be analyzed. 
  
  Default: \code{1}. }
  %%%%%%%%%%%%%%%%%%%%%%
  \item{exponent}{Exponent applied to averages before plotting; 
  the exponent should be positive.
  
  Default: \code{1}.
  }
  %%%%%%%%%%%%%%%%%%%%%%
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
  \item{minimum.level}{Minimum plot level of wavelet power averages considered. 
  
  Default: NULL (referring to minimum level observed).}
  \item{maximum.level}{Maximum plot level of wavelet power averages considered. 
  
  Default: NULL (referring to maximum level observed).}
  %%%%%%%%%%%%%%%%%%%%%%
  \item{label.avg.axis}{Label the axis of averages? Logical. 
  
  Default: \code{TRUE}.}
  \item{averagelab}{Label for the axis of averages. 
  
  Default: \code{"average wavelet power"}.}
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
  
  Default: \code{NULL}.}
  \item{lwd}{width of line of averages. 
  
  Default: \code{1}.}
  \item{col}{color of line of averages. 
  
  Default: \code{"black"}.}
  \item{lwd.axis}{line width of axes. 
  
  Default: \code{1}.}
  \item{verbose}{Print verbose output on the screen? Logical. 
  
  Default: \code{FALSE}.}
}
% \details{ }
% \value{ }
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

Tian, H., and Cazelles, B., 2012. \code{WaveletCo}. 
Available at \url{https://cran.r-project.org/src/contrib/Archive/WaveletCo/}, archived April 2013; accessed July 26, 2013.

Torrence C., and Compo G.P., 1998. 
A practical guide to wavelet analysis. 
Bulletin of the American Meteorological Society 79 (1), 61--78.
}
\author{Angi Roesch and Harald Schmidbauer; credits are also due to Huidong Tian and Bernard Cazelles}
% \note{ }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ 
\code{\link{analyze.wavelet}}, \code{\link{wt.image}}, \code{\link{wt.sel.phases}}, \code{\link{wt.phase.image}}, \code{\link{reconstruct}} 
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
        
my.data <- data.frame(x = x)

## Computation of wavelet power:
## a natural choice of 'dt' in the case of hourly data is 'dt = 1/24',
## resulting in one time unit equaling one day. 
## This is also the time unit in which periods are measured. 
my.wt <- analyze.wavelet(my.data, "x", loess.span = 0, 
                         dt = 1/24, dj = 1/20, 
                         lowerPeriod = 1/4,
                         make.pval = TRUE, n.sim = 10)

## Plot of wavelet power spectrum (with equidistant color breakpoints): 
wt.image(my.wt, color.key = "i", main = "wavelet power spectrum",
   legend.params = list(lab = "wavelet power levels (equidistant levels)"),
   periodlab = "period (days)")
## Note:
## The default time axis shows an index of given points in time, 
## which is the count of hours in our example.
         
## With time elapsed in days 
## (starting from 0 and proceeding in steps of 50 days) 
## instead of the (default) time index:
index.ticks  <- seq(1, series.length, by = 50*24)
index.labels <- (index.ticks-1)/24

## Insert your specification of time axis: 
wt.image(my.wt, color.key = "i", main = "wavelet power spectrum",
   legend.params = list(lab = "wavelet power levels (equidistant levels)"),
   periodlab = "period (days)", timelab = "time elapsed (days)",
   spec.time.axis = list(at = index.ticks, labels = index.labels))
         
## Plot of average wavelet power:
wt.avg(my.wt, siglvl = 0.05, sigcol = "red", periodlab = "period (days)")

## The same plot, but with enhanced symbol size, user-defined period axis, 
## and horizontal grid:
wt.avg(my.wt, siglvl = 0.05, sigcol = "red", sigcex = 1.3, 
   periodlab = "period (days)",
   spec.period.axis = list(at = c(1,8,32,128)), 
   periodtck = 1, periodtcl = NULL, 
   lwd = 1.5, lwd.axis = 0.25)
       
## Another style of the plot: 
## With user-defined period axis and axis of averages, 
## minimum and maximum plot levels of averages:

op <- par(no.readonly = TRUE)
par(cex.lab = 1.3, cex.axis = 1.1)
wt.avg(my.wt, siglvl = 0.05, sigcol = "red", sigcex = 1.3, 
   minimum.level = 0, maximum.level = 11,
   periodlab = "period (days)", 
   spec.period.axis = list(at = c(1,8,32,128)), 
   spec.avg.axis = list(at = 0:10),
   lwd = 1.5)
par(op)   
## Note:
## 'cex.axis' in 'par' controls for the size of axis tick labels, 
## while 'cex.lab' controls for the size of axis and legend labels. 
## Scaling by 'cex' would also affect 'sigcex'. 


## Please see also the examples in our guide booklet,
## URL http://www.hs-stat.com/projects/WaveletComp/WaveletComp_guided_tour.pdf.

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }