\name{wt.sel.phases}
\alias{wt.sel.phases}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot phases for selected periodic components of a single time series
}
\description{
This function plots the phases for selected periodic components of
a time series, which are provided by an object of class
\code{analyze.wavelet}.

Periodic components can be selected by specification of a single
Fourier period or of a period band. In the latter case, and in the
default case (no specification at all), there is an option to
average the phases across periods. Other options: restriction to the 
cone of influence, restriction to an area of significance (with respect 
to wavelet power).

(The time axis can be altered to give dates, see e.g. \code{wt.image}.)
}
\usage{
wt.sel.phases(WT, sel.period = NULL, sel.lower = NULL, sel.upper = NULL, 
      only.coi = F,
      only.sig = T, siglvl=0.05, 
      show.avg.phase = F, phase.avg.col = "black", 
      label.time.axis = T, show.date = F, date.format = NULL, timelab = NULL, 
      label.phase.axis = T, phaselab = NULL, main = NULL, sub = NULL, verbose = F)
}
\arguments{
  \item{WT}{ an object of class \code{analyze.wavelet}. }
  \item{sel.period}{ a single number which determines the (closest available) Fourier period to be selected. Default: NULL. }
  \item{sel.lower}{ a lower number which determines the lower (closest available) Fourier period to be selected if \code{sel.period} is NULL. Default: NULL. }
  \item{sel.upper}{ an upper number which determines the upper (closest available) Fourier period to be selected if \code{sel.period} is NULL. Default: NULL. } 
  \item{only.coi}{Exclude borders influenced by edge effects, i.e. include the cone of influence only? Logical. Default: FALSE.}
  \item{only.sig}{Use wavelet power significance to decide about the inclusion of (parts of) the phases' series? Logical. Default: TRUE.}
  \item{siglvl}{level of wavelet power significance. Default: 0.05.} 
  \item{show.avg.phase}{ Average phases over selected periods? (Effective only if a band of periods is selected.) Logical. Default: FALSE. }
  \item{phase.avg.col}{ Color for the plot of phase averages. Default: "black". }
  \item{label.time.axis}{ Label the time axis? Logical. Default: TRUE. }
  \item{show.date}{Show calendar dates? (Effective only if dates are available as rownames or by variable \code{date} in the data frame which has been analyzed using  \code{analyze.wavelet}.) Logical. Default: FALSE.}
  \item{date.format}{the format of date given as a character string, e.g. "\%Y-\%m-\%d \%H:\%M:\%S", or equivalently "\%F \%T"; see \code{strptime} for a list of implemented date conversion specifications. (If not specified, as.Date will be applied.) Default: NULL.}
  \item{timelab}{ Time axis label. Default: "time". }
  \item{label.phase.axis}{ Label the phase axis? Logical. Default: TRUE. }
  \item{phaselab}{ Phase axis label. Default: "phase". }
  \item{main}{ an overall title for the plot. Default: NULL. }
  \item{sub}{ a subtitle for the plot. Default: NULL. In this case, the selected period range will be given in the subtitle.}
  \item{verbose}{Print verbose output on the screen? Logical. Default: FALSE.}
}
% \details{ }
\value{A list of class \code{"sel.phases"} with the following elements:
   \item{Period}{the selected period (or period band)}
   \item{Phase}{time series of (average) phases at the selected period (or period band)} 
   \item{only.coi}{Is the influence of edge effects excluded? I.e. is the cone of influence used only?}  
   \item{only.sig}{Was wavelet power significance used in selection of phases?} 
   \item{siglvl}{level of wavelet power significance}
   \item{date}{time series of dates (if available)}
   \item{time.axis}{tick levels corresponding to the time steps used for wavelet transformation}
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
\code{\link{analyze.wavelet}}, \code{\link{wt.image}}, \code{\link{wt.avg}}, \code{\link{wt.phase.image}}, \code{\link{reconstruct}} 
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
     
my.date = seq(as.POSIXct("2014-10-14 00:00:00","\%F \%T"), by="hour", 
              length.out=series.length)     
my.data = data.frame(date=my.date, x=x)

my.w = analyze.wavelet(my.data, "x", loess.span=0, dt=1/24, dj=1/20, 
                       lowerPeriod=1/4, make.pval=T, n.sim=10)

## Plot of wavelet power spectrum (with equidistant color breakpoints):  
wt.image(my.w, color.key="i", timelab="time (days)", 
   legend.params=list(lab="wavelet power levels (equidistant levels)"))

## Select period 16 and plot corresponding phases across time:
wt.sel.phases(my.w, timelab="time (days)", sel.period=8)

## The same plot, but with calendar axis:
wt.sel.phases(my.w, timelab="", sel.period=8, 
              show.date=T, date.format="\%F \%T")

## In the following, no periods are selected;
## the plot shows average phases instead of individual phases: 
wt.sel.phases(my.w, timelab="time (days)", show.avg.phase=T)

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}