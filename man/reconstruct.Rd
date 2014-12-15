\name{reconstruct}
\alias{reconstruct}

\title{
Reconstruction of a (detrended) time series from output provided by an
object of class "analyze.wavelet" or "analyze.coherency"
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

(If the object provided as input is of class \code{analyze.coherency}, 
then the number or name of the time series can be specified.)  

Optional: plot of wavelets used for reconstruction, plot of 
reconstructed series against original (detrended)
series.

Output includes the original (detrended) and the reconstructed
time series, along with reconstruction waves and parameters.
}

\usage{
reconstruct(WT, my.series = 1, lvl = 0, 
      only.coi = F, 
      only.sig = T, siglvl = 0.05, 
      only.ridge = F, 
      sel.period = NULL, sel.lower = NULL, sel.upper = NULL,  
      rescale = T,
      plot.waves = F, plot.rec = T, 
      lty = 1, lwd = 1, col = 1:2, ylim = NULL,
      show.legend = T, 
      legend.coords = "topleft", legend.horiz=F, legend.text = NULL,
      label.time.axis = T, show.date = F, date.format = NULL, timelab = NULL, 
      main.waves = NULL, main.rec = NULL, main = NULL, verbose = T)
}

\arguments{
  \item{WT}{an object of class \code{analyze.wavelet} or \code{analyze.coherency}}
  \item{my.series}{In case class(WT) = \code{analyze.coherency}: number (1 or 2) or name of the series to be analyzed. Default: 1. }
  \item{lvl}{minimum level of wavelet power applied for inclusion of reconstruction waves. Default: 0.}
  \item{only.coi}{Exclude borders influenced by edge effects, i.e. include the cone of influence only? Logical. Default: FALSE.}
  \item{only.sig}{Use wavelet power significance at all? Logical. Default: TRUE.}
  \item{siglvl}{level of wavelet power significance applied for inclusion of reconstruction waves. Default: 0.05.} 
  \item{only.ridge}{Select ridge only? Logical. Default: FALSE.}
  \item{sel.period}{a vector of numbers to select Fourier periods (or closest available periods) and corresponding wavelets for the reconstruction. Default: NULL.}
  \item{sel.lower}{a number to define a lower Fourier period (or the closest available) for the selection of a band of wavelets for the reconstruction (only effective if sel.period is NULL). Default: NULL.}
  \item{sel.upper}{a number to define an upper Fourier period (or the closest available) for the selection of a band of wavelets for the reconstruction (only effective if sel.period is NULL). Default: NULL.}
  \item{rescale}{Should the reconstructed series finally be rescaled to attain the original (detrended) series' mean and standard deviation? Logical. Default: TRUE.}
  \item{plot.waves}{Should reconstruction waves be plotted? Logical. Default: FALSE.}
  \item{plot.rec}{Should the reconstructed series (together with the original (detrended) series) be plotted? Logical. Default: TRUE.}
  \item{lty}{parameter for the plot of original vs. reconstructed series: line type, e.g. 1:2. Default: 1.}
  \item{lwd}{parameter for the plot of original vs. reconstructed series: line width, e.g. 1:2. Default: 1.}
  \item{col}{parameter for the plot of original vs. reconstructed series: color of lines. Default: 1:2.}
  \item{ylim}{numeric vector of length 2, giving the time series coordinate range. Default: NULL.}
  \item{show.legend}{Include legend into the plot of original vs. reconstructed series? Logical. Default: TRUE.}
  \item{legend.coords}{coordinates to position the legend (as in function \code{legend}). Default: "topleft".}
  \item{legend.horiz}{ Set the legend horizontally rather than vertically? Logical. Default: FALSE.}
  \item{legend.text}{legend text. Default: c("original (detrended)", "reconstructed").}
  \item{label.time.axis}{Label the time axis? Logical. Default: TRUE.}
  \item{show.date}{Show calendar dates? (Effective only if dates are available as rownames or by variable \code{date} in the data frame which has been analyzed using  \code{analyze.wavelet}.) Logical. Default: FALSE.}
  \item{date.format}{the format of date given as a character string, e.g. "\%Y-\%m-\%d", or equivalently "\%F"; see \code{strptime} for a list of implemented date conversion specifications. (If not specified, as.Date will be applied.) Default: NULL.}
  \item{timelab}{Time axis label. Default: "time".}
  \item{main.waves}{an overall title for the plot of reconstruction waves. Default: NULL.}
  \item{main.rec}{an overall title for the plot of original vs. reconstructed series. Default: NULL.}
  \item{main}{an overall title for both plots. Default: NULL.}
  \item{verbose}{Print verbose output on the screen? Logical. Default: TRUE.}
}
% \details{ }
\value{A list of class \code{reconstruct} with the following elements:
   \item{series}{a data frame with the following columns
  
            \tabular{rllll}{
                     \tab date             \tab : \tab the calendar date (if available as column in my.data) \cr
                     \tab <x>              \tab : \tab series <x>, with original name retained \cr
                     \tab                  \tab : \tab (detrended, if loess.span != 0) \cr
                     \tab <x>.trend        \tab : \tab the trend series (if loess.span != 0) \cr
                     \tab <x>.r            \tab : \tab the reconstructed (detrended) series
                     } 
                     
             Row names are taken over from WT.  
   
   }
   \item{rec.waves}{data frame of scaled waves used for reconstruction}
   \item{loess.span}{parameter alpha in loess controlling the degree of time series smoothing, if the time series was detrended; no detrending if loess.span=0.}
   \item{lvl}{level which the wavelet power should have at least for waves (wave segments) to be included in the reconstruction. Default: 0.}
   \item{only.coi}{Is the influence of edge effects excluded? I.e. is the cone of influence used only?}  
   \item{only.sig}{Was wavelet power significance used in reconstruction?} 
   \item{siglvl}{level of wavelet power significance}
   \item{only.ridge}{Select ridge only? Logical.}  
   \item{rnum.used}{the vector of Fourier period numbers used for reconstruction }  
   \item{rescale}{Was the reconstructed series rescaled according to the mean and standard deviation taken from the original (detrended) series?}
   \item{dt}{time resolution, i.e. sampling resolution on time domain, 1/dt = number of intervals per time step.}
   \item{dj}{frequency resolution, i.e. sampling resolution on frequency domain, 1/dj = number of suboctaves (voices per octave).}
   \item{Period}{the Fourier periods (in time units)}
   \item{Scale}{the scales}
   \item{nc}{number of columns/time steps}
   \item{nr}{number of rows/scales}
   \item{axis.1}{tick levels corresponding to time steps}
   \item{axis.2}{tick levels corresponding to Fourier periods = log2(Period)}
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
\code{\link{analyze.wavelet}}, \code{\link{wt.avg}}, \code{\link{wt.image}}, \code{\link{wt.sel.phases}}, \code{\link{wt.phase.image}} 
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
 xlab="time (days)", 
 ylab="hourly data", main="a series of hourly data with periods of 1, 8, 32, and 128 days")
  
my.data = data.frame(x=x)

my.w = analyze.wavelet(my.data, "x", 
                       loess.span=0, 
                       dt=1/24, dj=1/20, 
                       lowerPeriod = 1/4, 
                       make.pval=T, n.sim=10)

## Plot of wavelet power spectrum (with equidistant color breakpoints):  
wt.image(my.w, color.key="interval", legend.params=list(lab="wavelet power levels"))

## Reconstruction of the time series, including significant components only: 
reconstruct(my.w, timelab="time (days)")
## The same reconstruction, but showing wave components first:
reconstruct(my.w, timelab="time (days)", plot.waves=T)

## Reconstruction, including all components whether significant or not: 
reconstruct(my.w, timelab="time (days)", only.sig=F)

## Reconstruction, including significant components, but selected periods only:  
reconstruct(my.w, timelab="time (days)", sel.period=c(1,8,32,128))

## Reconstruction, including significant components, but the ridge only:
reconstruct(my.w, timelab="time (days)", only.ridge=T)

## See the periods involved:
my.rec = reconstruct(my.w, timelab="time (days)", only.ridge=T)
print(my.rec$Period[my.rec$rnum.used])
## The original and reconstructed time series can be retrieved as well:
plot(my.rec$series$x, type="l", xlab="time (days)", ylab="")
lines(my.rec$series$x.r, col="red")
legend("topleft", legend=c("original","reconstructed"), lty=1, col=c("black","red"))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }