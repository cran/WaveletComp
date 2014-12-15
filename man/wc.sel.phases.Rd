\name{wc.sel.phases}
\alias{wc.sel.phases}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Comparison plot of phases for selected periodic components of two time
series
}
\description{
This function plots the phases for periodic components of two time series, which are provided by an object
of class \code{analyze.coherency}.

Periodic components can be selected by specification of a single
Fourier period or of a period band. In the latter case, and in the
default case (no specification at all), phases are averaged across
periods for each time series. Other options: restriction to the 
cone of influence, restriction to an area of significance (with respect 
to cross-wavelet power, wavelet coherence or individual wavelet power). 
Phase differences (i.e. angles, smoothed or not
smoothed) can be added to the plot.

(The time axis can be altered to display dates, see e.g. \code{wt.image}.)
}
\usage{
wc.sel.phases(WC, sel.period = NULL, sel.lower = NULL, sel.upper = NULL, 
      only.coi = F, 
      only.sig = T, which.sig = "wp", siglvl=0.05, 
      phase.cols = c("red", "blue"), 
      show.Angle = T, use.sAngle = F, Angle.col = "black", 
      show.legend = T, legend.coords = "topleft", legend.horiz = T,
      label.time.axis = T, show.date = F, date.format = NULL, timelab = NULL, 
      label.phase.axis = T, phaselab = NULL, 
      phaselim = c(-pi,pi+show.legend*ifelse(legend.horiz,0.8,2)),      
      main = NULL, sub = NULL,
      verbose = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{WC}{ an object of class \code{analyze.wavelet}. }
  \item{sel.period}{ a single number which determines the (closest available) Fourier period to be selected. Default: NULL. }
  \item{sel.lower}{ a lower number which determines the lower (closest available) Fourier period to be selected if \code{sel.period} is NULL. Default: NULL. }
  \item{sel.upper}{ an upper number which determines the upper (closest available) Fourier period to be selected if \code{sel.period} is NULL. Default: NULL. }
  \item{only.coi}{Exclude borders influenced by edge effects, i.e. include the cone of influence only? Logical. Default: FALSE.}
  \item{only.sig}{Use cross-wavelet power or coherence significance to decide about the inclusion of (parts of) the phases' series? Logical. Default: TRUE.}
  \item{which.sig}{Which spectrum should significance refer to?
                     
                       \tabular{rlll}{
                                \tab "wp" \tab : \tab cross-wavelet power (default) \cr
                                \tab "wc" \tab : \tab wavelet coherence \cr
                                \tab "wt" \tab : \tab individual wavelet power
                                }
                                
                        Default: "wp"
 
  }
  \item{siglvl}{level of significance. Default: 0.05.} 
  \item{phase.cols}{ a vector of two colors for the plot of (average) phases referring to the two time series. Default: c("red","blue"). }
  \item{show.Angle}{ Show the (average) phase difference (the Angle) between the two series? Logical. Default: TRUE. }
  \item{use.sAngle}{ Use smoothed version of phase difference? Logical. Default: FALSE. }
  \item{Angle.col}{ Color for the plot of Angles. Default: "black". }
  \item{show.legend}{ Include legend? Logical. Default: TRUE. }
  \item{legend.coords}{ Coordinates to position the legend (with the same options as given in function \code{legend}). Default: "topleft". }
  \item{legend.horiz}{ Set the legend horizontally rather than vertically? Logical. Default: TRUE.}
  \item{label.time.axis}{ Label the time axis? Logical. Default: TRUE. }
  \item{show.date}{Show calendar dates? (Effective only if dates are available as rownames or as variable \code{date} in the data frame analyzed using \code{analyze.coherency}.) Logical. Default: FALSE.}
  \item{date.format}{the format of date given as a character string, e.g. "\%Y-\%m-\%d", or equivalently "\%F"; see \code{strptime} for a list of implemented date conversion specifications. (If not specified, as.Date will be applied.) Default: NULL.}
  \item{timelab}{ Time axis label. Default: "time". }
  \item{label.phase.axis}{Label the phase axis? Logical. Default: TRUE. }
  \item{phaselab}{ Phase axis label. Default: "phase". }   
  \item{phaselim}{numeric vector of length 2, giving the phase coordinate range.
  Default: c(-pi,pi+0.8) (+0.8 in order to accomodate the horizontal legend, +2 in case of a vertical legend).}
  \item{main}{ an overall title for the plot. Default: NULL.}
  \item{sub}{ a subtitle for the plot. Default: NULL. In this case, the selected period range will be given in the subtitle.}
  \item{verbose}{Print verbose output on the screen? Logical. Default: FALSE.}
}
% \details{ }
\value{A list of class \code{sel.phases} with the following elements:
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
   \item{date}{time series of dates (if available)}
   \item{time.axis}{tick levels corresponding to the time steps used for wavelet transformation}
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
\code{\link{analyze.wavelet}}, \code{\link{wt.image}}, \code{\link{wt.avg}}, \code{\link{wt.phase.image}}, \code{\link{reconstruct}}
}
\examples{
\dontrun{
## The following example is adopted from Veleda et al, 2012

add.noise=TRUE

series.length = 3*128*24
x1 = periodic.series(start.period = 1*24, length = series.length)
x2 = periodic.series(start.period = 2*24, length = series.length)
x3 = periodic.series(start.period = 4*24, length = series.length)
x4 = periodic.series(start.period = 8*24, length = series.length)
x5 = periodic.series(start.period = 16*24, length = series.length)
x6 = periodic.series(start.period = 32*24, length = series.length)
x7 = periodic.series(start.period = 64*24, length = series.length)
x8 = periodic.series(start.period = 128*24, length = series.length)

x = x1 + x2 + x3 + x4 + 3*x5 + x6 + x7 + x8 
y = x1 + x2 + x3 + x4 + 3*x5 + x6 + 3*x7 + x8

if (add.noise == TRUE){
    x = x + rnorm(length(x))
    y = y + rnorm(length(y))
}

my.data = data.frame(x=x, y=y)
ts.plot(ts(my.data$x, start=0, frequency=24), 
     ts(my.data$y, start=0, frequency=24), type="l", col=1:2, 
     xlab="time (days)", ylab="hourly data", 
     main="a series of hourly data with periods of 1, 2, 4, 8, 16, 32, 64, and 128 days",
     sub="(different amplitudes at periods 16 and 64)")
legend("topright", legend=c("x","y"), col=1:2, lty=1)

## computation of cross-wavelet power and wavelet coherency
my.wc = analyze.coherency(my.data, c("x","y"), loess.span=0, 
                           dt=1/24, dj=1/20,
                           window.size.t=1, window.size.s=1/2, 
                           lowerPeriod=1/4, 
                           make.pval=T, n.sim=10)

## plot of cross-wavelet power
wc.image(my.wc, timelab="time (days)", periodlab="period (days)",
         main="cross-wavelet power")

## Select period 64 and compare plots of corresponding phases, including the 
## phase differences (angles) in their non-smoothed (default) version:
wc.sel.phases(my.wc, timelab="time (days)", sel.period=64, show.Angle=T)

## In the following, no periods are selected. In this case, instead of individual phases 
## the plot shows average phases: 
wc.sel.phases(my.wc, timelab="time (days)")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }