\name{analyze.coherency}
\alias{analyze.coherency}

\title{
Computation of the cross-wavelet power and wavelet coherence spectrum
of two time series
}

\description{
The two time series are selected from an input data frame by specifying either their names 
or their column numbers. Optionally, the time series are detrended, using \code{loess} with parameter
\code{loess.span}. Internally, the series will be standardized before they undergo wavelet transformation.

The cross-wavelet power spectrum is computed applying the Morlet wavelet. 
P-values to test the null hypothesis that a period (within \code{lowerPeriod} and \code{upperPeriod}) 
is irrelevant at a certain time are calculated if desired; this is accomplished with the help of a
simulation algorithm. There is a selection of models from which to choose the alternative hypothesis.
The selected model will be fitted to the data and simulated according to estimated parameters
in order to provide surrogate time series.

For the computation of wavelet coherence, a variety of filtering methods is provided, 
with flexible window parameters.

Wavelet transformation, as well as p-value computations, are carried out by calling subroutine \code{wc}.

The name and parts of the layout of subroutine \code{wc} were inspired by a similar function 
developed by Huidong Tian and Bernard Cazelles (archived R package \code{WaveletCo}). 
The basic concept of the simulation algorithm and of ridge determination build on ideas 
developed by these authors. The major part of the code for the computation of the cone of influence
and the code for Fourier-randomized surrogate time series has been adopted from Huidong Tian. 
The implementation of a choice of filtering windows for the computation of the wavelet coherence 
was inspired by Luis Aguiar-Conraria and Maria Joana Soares (\code{GWPackage}).

Cross-wavelet and coherence computation, the simulation algorithm and ridge determination build heavily 
on the use of matrices in order to minimize computation time in R. 

This function provides a broad variety of final as well as intermediate results 
which can be further analyzed in detail. 
}
\usage{
analyze.coherency(my.data, my.pair = c(1, 2), loess.span = 0.75, 
                  dt = 1, dj = 1/20, 
                  lowerPeriod = 2*dt, 
                  upperPeriod = floor(nrow(my.data)/3)*dt, 
                  window.type.t = 1, window.type.s = 1, 
                  window.size.t = 5, window.size.s = 1/4, 
                  make.pval = TRUE, method = "white.noise", params = NULL, 
                  n.sim = 100, 
                  date.format = NULL, date.tz = NULL, 
                  verbose = TRUE)
}

\arguments{
  \item{my.data}{data frame of time series (including header, and dates as row names or as separate column 
  named \code{"date"} if available)}
  \item{my.pair}{pair of names or column indices indicating the series to be analyzed, 
  e.g. \code{c(1,2)}, \code{c(2,1)}, \code{c("dji","ftse")}. 
  
  Default: \code{c(1,2)}.
  }
  \item{loess.span}{parameter \code{alpha} in \code{loess} controlling the degree of time series smoothing, 
  if the time series are to be detrended; no detrending if \code{loess.span = 0}. 
  
  Default: \code{0.75}.
  }
  %%%%%%%%%%%%%%%%%
  \item{dt}{time resolution, i.e. sampling resolution in the time domain, \code{1/dt} = number of observations per time unit. 
  For example: a natural choice of \code{dt} in case of hourly data is \code{dt = 1/24}, 
  resulting in one time unit equaling one day. This is also the time unit in which periods are measured.
  If \code{dt = 1}, the time interval between two consecutive observations will equal one time unit.
  
  Default: \code{1}.
  }
  \item{dj}{frequency resolution, i.e. sampling resolution in the frequency domain, \code{1/dj} = number of 
  suboctaves (voices per octave).
  
  Default: \code{1/20}.
  }
  %%%%%%%%%%%%%%%%%
  \item{lowerPeriod}{lower Fourier period (measured in time units determined by \code{dt}, 
  see the explanations concerning \code{dt}) for wavelet decomposition.\cr 
  If \code{dt = 1}, the minimum admissible value is 2. 
  
  Default: \code{2*dt}.}
  \item{upperPeriod}{upper Fourier period (measured in time units determined by \code{dt}, 
  see the explanations concerning \code{dt}) for wavelet decomposition. 
  
  Default: \code{(floor of one third of time series length)*dt}.
  } 
  %%%%%%%%%%%%%%%%%
  \item{window.type.t}{type of window for smoothing in time direction; select from:  
                 \tabular{rllll}{
                     \tab  \code{0} \tab (\code{"none"}) \tab : \tab no smoothing in time direction \cr
                     \tab  \code{1} \tab (\code{"bar"})  \tab : \tab Bartlett \cr
                     \tab  \code{2} \tab (\code{"tri"})  \tab : \tab Triangular (Non-Bartlett) \cr
                     \tab  \code{3} \tab (\code{"box"})  \tab : \tab Boxcar (Rectangular, Dirichlet)\cr
                     \tab  \code{4} \tab (\code{"han"})  \tab : \tab Hanning \cr
                     \tab  \code{5} \tab (\code{"ham"})  \tab : \tab Hamming \cr
                     \tab  \code{6} \tab (\code{"bla"})  \tab : \tab Blackman 
                     } 
                 Default: \code{1} = \code{"bar"}.
                 }
  \item{window.type.s}{type of window for smoothing in scale (period) direction; select from:
                 \tabular{rllll}{
                     \tab  \code{0} \tab (\code{"none"}) \tab : \tab no smoothing in scale (period) \cr 
                     \tab           \tab                 \tab   \tab direction \cr
                     \tab  \code{1} \tab (\code{"bar"})  \tab : \tab Bartlett \cr
                     \tab  \code{2} \tab (\code{"tri"})  \tab : \tab Triangular (Non-Bartlett) \cr
                     \tab  \code{3} \tab (\code{"box"})  \tab : \tab Boxcar (Rectangular, Dirichlet)\cr
                     \tab  \code{4} \tab (\code{"han"})  \tab : \tab Hanning \cr
                     \tab  \code{5} \tab (\code{"ham"})  \tab : \tab Hamming \cr
                     \tab  \code{6} \tab (\code{"bla"})  \tab : \tab Blackman 
                     } 
                 Default: \code{1} = \code{"bar"}.
  }
  \item{window.size.t}{size of the window used for smoothing in time direction, measured in time units 
  determined by \code{dt}, see the explanations concerning \code{dt}. 
  
  Default: \code{5}, which together with \code{dt = 1} defines a window of length \code{5*(1/dt)} = 5, 
  equaling 5 observations (observation epochs). Windows of even-numbered sizes are extended by 1.
  }
  \item{window.size.s}{size of the window used for smoothing in scale (period) direction in units of \code{1/dj}. 
  
  Default: \code{1/4}, which together with \code{dj = 1/20} defines a window of length \code{(1/4)*(1/dj)} = 5. 
  Windows of even-numbered sizes are extended by 1.
  } 
  %%%%%%%%%%%%%%%%%
  \item{make.pval}{Compute p-values? Logical.
  
  Default: \code{TRUE}.}
  %%%%%%%%%%%%%%%%% 
  \item{method}{the method of generating surrogate time series; select from:
  
                \tabular{rlll}{
                    \tab \code{"white.noise"}  \tab : \tab white noise \cr
                    \tab \code{"shuffle"}      \tab : \tab shuffling the given time series \cr
                    \tab \code{"Fourier.rand"} \tab : \tab time series with a similar spectrum \cr
                    \tab \code{"AR"}           \tab : \tab AR(p) \cr  
                    \tab \code{"ARIMA"}        \tab : \tab ARIMA(p,0,q) 
                    }
  Default: \code{"white.noise"}.
  }
  \item{params}{a list of assignments between methods (AR, and ARIMA) and lists of parameter values 
  applying to surrogates. Default: \code{NULL}. 
  
  Default includes two lists named \code{AR} and \code{ARIMA}:
  \itemize{  
  \item \code{AR = list(...)}, a list containing one single element:
        \tabular{rlll}{
        \tab \code{p} \tab : \tab AR order.\cr
        \tab          \tab   \tab Default: \code{1}.
                }
  \item \code{ARIMA = list(...)}, a list of six elements:            
        \tabular{rlll}{
        \tab \code{p}            \tab : \tab  AR order. \cr
        \tab                     \tab   \tab  Default: \code{1}.\cr
        \tab \code{q}            \tab : \tab  MA order. \cr
        \tab                     \tab   \tab  Default: \code{1}.\cr
        \tab \code{include.mean} \tab : \tab  Include a mean/intercept term? \cr
        \tab                     \tab   \tab  Default: \code{TRUE}.\cr
        \tab \code{sd.fac}       \tab : \tab  magnification factor to boost the \cr 
        \tab                     \tab   \tab  residual standard deviation. \cr  
        \tab                     \tab   \tab  Default: \code{1}.\cr
        \tab \code{trim}         \tab : \tab  Simulate trimmed data? \cr
        \tab                     \tab   \tab  Default: \code{FALSE}.\cr
        \tab \code{trim.prop}    \tab : \tab  high/low trimming proportion. \cr
        \tab                     \tab   \tab  Default: \code{0.01}.\cr
                }}
  }      
  %%%%%%%%%%%%%%%%%
  \item{n.sim}{number of simulations. 
  
  Default: \code{100}.}
  %%%%%%%%%%%%%%%%%
  \item{date.format}{optional, and for later reference: the format of calendar date 
  (if available in the input data frame) given as a character string, e.g. \code{"\%Y-\%m-\%d"}, 
  or equivalently \code{"\%F"}; see \code{strptime} for a list of implemented date conversion specifications.
  Explicit information given here will be overwritten by any later specification given in 
  e.g. \code{wc.image}. 
  If unspecified, date formatting will be attempted according to \code{as.Date}. 
  
  Default: \code{NULL}.
  }
  \item{date.tz}{optional, and for later reference: a character string specifying the time zone of calendar date 
  (if available in the input data frame); see \code{strptime}. 
  Explicit information given here will be overwritten by any specification given in 
  e.g. \code{wc.image}. 
  If unspecified, \code{""} (the local time zone) will be used. 
  
  Default: \code{NULL}.
  }
  %%%%%%%%%%%%%%%%%  
  \item{verbose}{Print verbose output on the screen? Logical. 
  
  Default: \code{TRUE}.}
}
% \details{ }
\value{A list of class \code{"analyze.coherency"} with elements of different dimensions.
%%%%%%%%%%%%%%%%%
The elements of matrix type, namely:
\itemize{
\item \code{Wave.xy}, \code{Angle}, \code{sWave.xy}, \code{sAngle}, 
\item \code{Power.xy}, \code{Power.xy.pval}, 
\item \code{Coherency}, \code{Coherence}, \code{Coherence.pval}, 
\item \code{Wave.x}, \code{Wave.y}, \code{Phase.x}, \code{Phase.y}, \code{Ampl.x}, \code{Ampl.y}, 
\item \code{Power.x}, \code{Power.y}, \code{Power.x.pval}, \code{Power.y.pval}, \code{sPower.x}, \code{sPower.y}, 
\item \code{Ridge.xy}, \code{Ridge.co}, \code{Ridge.x}, \code{Ridge.y},} 
have the following structure:\cr 
columns correspond to observations (observation epochs; "epoch" meaning point in time), 
rows correspond to scales (Fourier periods) whose values are given in \code{Scale} (\code{Period}).
%%%%%%%%%%%%%%%%%
Here is a detailed list of all elements:
  \item{series}{a data frame with the following columns:
  
            \tabular{rlll}{
                     \tab date                 \tab : \tab the calendar date \cr
                     \tab                      \tab   \tab (if available as column in \code{my.data}) \cr
                     \tab <x>, <y>             \tab : \tab the two series which have been analyzed \cr
                     \tab                      \tab   \tab (detrended, if \code{loess.span != 0}; \cr
                     \tab                      \tab   \tab original names retained) \cr
                     \tab <x>.trend, <y>.trend \tab : \tab the two trend series \cr
                     \tab                      \tab   \tab (included if \code{loess.span != 0}) 
                     } 
                     
  Row names are taken over from \code{my.data}, and so are dates if given as row names.        
  }    
  \item{loess.span}{parameter \code{alpha} in \code{loess} controlling the degree of time series smoothing 
  if the time series were detrended; no detrending if \code{loess.span = 0}}
  %%%%%%%%%%%%%%%%%
  \item{dt}{time resolution, i.e. sampling resolution in the time domain, \code{1/dt} = number of observations per time unit}
  \item{dj}{frequency resolution, i.e. sampling resolution in the frequency domain, \code{1/dj} = number of suboctaves (voices per octave)}
  %%%%%%%%%%%%%%%%%
  \item{Wave.xy}{(complex-valued) cross-wavelet transform (analogous to Fourier cross-frequency spectrum, and to the covariance in statistics)} 
  \item{Angle}{phase difference, i.e. phase lead of <x> over <y> (= \code{phase.x} minus \code{phase.y})}
  \item{sWave.xy}{smoothed (complex-valued) cross-wavelet transform}
  \item{sAngle}{phase difference, i.e. phase lead of <x> over <y>, affected by smoothing} 
  %%%%%%%%%%%%%%%%%
  \item{Power.xy}{cross-wavelet power (analogous to Fourier cross-frequency power spectrum)}
  \item{Power.xy.avg}{average cross-wavelet power in the frequency domain (averages over time)} 
  \item{Power.xy.pval}{p-values of cross-wavelet power}
  \item{Power.xy.avg.pval}{p-values of average cross-wavelet power}
  %%%%%%%%%%%%%%%%%
  \item{Coherency}{the (complex-valued) wavelet coherency of series <x> over series <y> in the time/frequency domain, affected by smoothing (analogous to Fourier coherency, and to the coefficient of correlation in statistics)}       
  \item{Coherence}{wavelet coherence (analogous to Fourier coherence, and to the coefficient of determination in statistics (affected by smoothing)}
  \item{Coherence.avg}{average wavelet coherence in the frequency domain (averages across time)}
  \item{Coherence.pval}{p-values of wavelet coherence}
  \item{Coherence.avg.pval}{p-values of average wavelet coherence}
  %%%%%%%%%%%%%%%%%
  \item{Wave.x, Wave.y}{(complex-valued) wavelet transforms of series <x> and <y>}
  \item{Phase.x, Phase.y}{phases of series <x> and <y>}
  \item{Ampl.x, Ampl.y}{amplitudes of series <x> and <y>}
  \item{Power.x, Power.y}{wavelet power of series <x> and <y>}
  \item{Power.x.avg, Power.y.avg}{average wavelet power of series <x> and <y>, averages across time} 
  \item{Power.x.pval, Power.y.pval}{p-values of wavelet power of series <x> and <y>}
  \item{Power.x.avg.pval, Power.y.avg.pval}{p-values of average wavelet power of series <x> and <y>}
  \item{sPower.x, sPower.y}{smoothed wavelet power of series <x> and <y>} 
  %%%%%%%%%%%%%%%%%
  \item{Ridge.xy}{ridge of cross-wavelet power, in the form of a matrix of \code{0}s and \code{1}s}
  \item{Ridge.co}{ridge of wavelet coherence}
  \item{Ridge.x, Ridge.y}{power ridges of series <x> and <y>}
  %%%%%%%%%%%%%%%%%
  \item{Period}{the Fourier periods 
  (measured in time units determined by \code{dt}, see the explanations concerning \code{dt})}
  \item{Scale}{the scales (the Fourier periods divided by the Fourier factor)}
  %%%%%%%%%%%%%%%%%
  \item{nc}{number of columns = number of observations = number of observation epochs; "epoch" meaning point in time}
  \item{nr}{number of rows = number of scales (Fourier periods)}
  %%%%%%%%%%%%%%%%%
  \item{coi.1, coi.2}{borders of the region where the wavelet transforms are not influenced by edge effects (cone of influence). 
  The coordinates of the borders are expressed in terms of internal axes \code{axis.1} and \code{axis.2}.}
  %%%%%%%%%%%%%%%%%
  \item{axis.1}{tick levels corresponding to the time steps used for (cross-)wavelet transformation: \code{1, 1+dt, 1+2dt, ...}.
  The default time axis in plot functions provided by \code{WaveletComp} is determined by observation epochs, however; "epoch" meaning point in time. }
  \item{axis.2}{tick levels corresponding to the log of Fourier periods: \code{log2(Period)}. This determines the period axis in plot functions provided by \code{WaveletComp}.}
  %%%%%%%%%%%%%%%%%
  \item{date.format}{the format of calendar date (if available)}
  \item{date.tz}{the time zone of calendar date (if available)}
}
\references{ 

Aguiar-Conraria L., and Soares M.J., 2011. 
Business cycle synchronization and the Euro: A wavelet analysis. 
Journal of Macroeconomics 33 (3), 477--489.

Aguiar-Conraria L., and Soares M.J., 2011. 
The Continuous Wavelet Transform: A Primer. 
NIPE Working Paper Series 16/2011.

Aguiar-Conraria L., and Soares M.J., 2012. \code{GWPackage}. 
Available at \url{https://sites.google.com/site/aguiarconraria/joanasoares-wavelets}; accessed September 4, 2013.

Cazelles B., Chavez M., Berteaux, D., Menard F., Vik J.O., Jenouvrier S., and Stenseth N.C., 2008.
Wavelet analysis of ecological time series.
Oecologia 156, 287--304.

Liu P.C., 1994. 
Wavelet spectrum analysis and ocean wind waves. 
In: Foufoula-Georgiou E., and Kumar P., (eds.), Wavelets in Geophysics, Academic Press, San Diego, 151--166.

Tian, H., and Cazelles, B., 2012. \code{WaveletCo}. 
Available at \url{https://cran.r-project.org/src/contrib/Archive/WaveletCo/}, archived April 2013; accessed July 26, 2013.

Torrence C., and Compo G.P., 1998. 
A practical guide to wavelet analysis.
Bulletin of the American Meteorological Society 79 (1), 61--78.

Veleda D., Montagne R., and Araujo M., 2012. 
Cross-Wavelet Bias Corrected by Normalizing Scales. 
Journal of Atmospheric and Oceanic Technology 29, 1401--1408.
}
\author{Angi Roesch and Harald Schmidbauer; 
credits are also due to Huidong Tian, Bernard Cazelles, Luis Aguiar-Conraria, and Maria Joana Soares.}
% \note{ }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ 
\code{\link{wc.image}}, \code{\link{wc.avg}}, \code{\link{wc.sel.phases}}, \code{\link{wc.phasediff.image}}, 
\code{\link{wt.image}}, \code{\link{wt.avg}}, \code{\link{wt.sel.phases}}, \code{\link{wt.phase.image}}, \code{\link{reconstruct}} 
}
\examples{
\dontrun{
## The following example is modified from Veleda et al, 2012:

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
## There is an option to store the date format and time zone as additional
## parameters within object 'my.wc' for later reference. 

my.wc <- analyze.coherency(my.data, c("x","y"), 
                           loess.span = 0, 
                           dt = 1/24, dj = 1/20, 
                           window.size.t = 1, window.size.s = 1/2, 
                           lowerPeriod = 1/4,
                           make.pval = TRUE, n.sim = 10,
                           date.format = "\%F \%T", date.tz = "")
## Note:                           
## By default, Bartlett windows are used for smoothing in order to obtain
## the wavelet coherence spectrum; window lengths in this example:
## 1*24 + 1 = 25 observations in time direction,
## (1/2)*20 + 1 = 11 units in scale (period) direction.                             
                         
## Plot of cross-wavelet power 
## (with color breakpoints according to quantiles):
wc.image(my.wc, main = "cross-wavelet power spectrum, x over y",
   legend.params = list(lab = "cross-wavelet power levels"),
   periodlab = "period (days)")
   
## The same plot, now with calendar axis
## (according to date format stored in 'my.wc'):
wc.image(my.wc, main = "cross-wavelet power spectrum, x over y",
   legend.params = list(lab = "cross-wavelet power levels"),
   periodlab = "period (days)", show.date = TRUE)   
         
## Plot of average cross-wavelet power:
wc.avg(my.wc, siglvl = 0.05, sigcol = 'red', 
   periodlab = "period (days)")

## Plot of wavelet coherence 
## (with color breakpoints according to quantiles):
wc.image(my.wc, which.image = "wc",  main = "wavelet coherence, x over y", 
   legend.params = list(lab = "wavelet coherence levels", 
                        lab.line = 3.5, label.digits = 3),
   periodlab = "period (days)")
         
## plot of average coherence:
wc.avg(my.wc, which.avg = "wc", 
   siglvl = 0.05, sigcol = 'red', 
   legend.coords = "topleft", 
   periodlab = "period (days)")

## Please see our guide booklet for further examples:
## URL http://www.hs-stat.com/projects/WaveletComp/WaveletComp_guided_tour.pdf.

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }