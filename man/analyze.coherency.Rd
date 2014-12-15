\name{analyze.coherency}
\alias{analyze.coherency}

\title{
Computation of the cross-wavelet power and wavelet coherence spectrum
of two time series
}

\description{
The two time series are selected from an input data frame by
specifying either their names or their column numbers. Optionally,
the time series are detrended, using \code{loess} with parameter
\code{loess.span}. Internally, the series will be
standardized before they undergo wavelet transformation.

The cross-wavelet power spectrum is computed applying the Morlet
wavelet. P-values to test the null hypothesis that a period (within lowerPeriod and upperPeriod) is irrelevant at a certain time
are calculated if desired; this is accomplished with the help of a
simulation algorithm. There is a selection of models from which to
choose the alternative hypothesis.
The selected model will be fitted to
the data and simulated according to estimated parameters in order
to provide surrogate time series.

For the computation of wavelet coherence, a variety of filtering
methods is provided, with flexible window parameters.

Wavelet transformation, as well as p-value computations, are
carried out by calling subroutine \code{wc}.

The name and parts of the layout of subroutine \code{wc} were inspired by a similar function developed by 
Huidong Tian and Bernard Cazelles (archived R package \code{WaveletCo}). The basic concept of the simulation algorithm, and of ridge determination build on ideas developed by these authors. The major part of the code for the computation of the cone of influence, and the code for Fourier-randomized surrogate time series
has been adopted from Huidong Tian. 
The implementation of a choice of filtering windows for the computation of
the wavelet coherence was inspired by Luis Aguiar-Conraria and Maria Joana Soares (\code{GWPackage}).

Cross-wavelet and coherence computation, the simulation algorithm, and ridge determination build heavily on the use of matrices in order to minimize computation time in R. 

This function provides a broad variety of final as well as intermediate results 
which can be further analyzed in detail. 
}
\usage{
analyze.coherency(my.data, my.pair = c(1, 2), loess.span = 0.75, dt = 1, dj = 1/20, 
                  lowerPeriod = 2*dt, upperPeriod = floor(nrow(my.data)/3)*dt, 
                  window.type.t = 1, window.type.s = 1, 
                  window.size.t = 5, window.size.s = 1/4, 
                  make.pval = T, method = "white.noise", params = NULL, 
                  n.sim = 100, verbose = T)
}

\arguments{
  \item{my.data}{data frame of time series (including header, and dates as row names or as separate column named "date" if available)}
  \item{my.pair}{pair of names or column indices indicating series x and y to be analyzed, e.g. c(1,2), c(2,1), c("dji","ftse"). 
  
  Default: c(1,2).
  }
  \item{loess.span}{parameter alpha in loess controlling the degree of time series smoothing, if the time series is to be detrended; no detrending if loess.span=0. 
  
  Default: 0.75.
  }
  \item{dt}{time resolution, i.e. sampling resolution on time domain, 1/dt = number of intervals per time step. Default: 1.}
  \item{dj}{frequency resolution, i.e. sampling resolution on frequency domain, 1/dj = number of suboctaves (voices per octave). Default: 1/20.}
  \item{lowerPeriod}{lower Fourier period (in time units) for wavelet decomposition. Default: 2*dt.}
  \item{upperPeriod}{upper Fourier period (in time units) for wavelet decomposition. 
  
                     Default: (floor of one third of time series length)*dt.} 
  \item{window.type.t}{type of window for smoothing in time direction, select from:  
                 \tabular{rllll}{
                     \tab  0 \tab ("none") \tab : \tab no smoothing in time direction \cr
                     \tab  1 \tab ("bar")  \tab : \tab Bartlett \cr
                     \tab  2 \tab ("tri")  \tab : \tab Triangular (Non-Bartlett) \cr
                     \tab  3 \tab ("box")  \tab : \tab Boxcar (Rectangular, Dirichlet)\cr
                     \tab  4 \tab ("han")  \tab : \tab Hanning \cr
                     \tab  5 \tab ("ham")  \tab : \tab Hamming \cr
                     \tab  6 \tab ("bla")  \tab : \tab Blackman 
                     } 
                 Default: 1 = "bar".
                 }
  \item{window.type.s}{type of window for smoothing in scale (period) direction, select from:
                 \tabular{rllll}{
                     \tab  0 \tab ("none") \tab : \tab no smoothing in scale (period) direction \cr
                     \tab  1 \tab ("bar") \tab : \tab Bartlett \cr
                     \tab  2 \tab ("tri") \tab : \tab Triangular (Non-Bartlett) \cr
                     \tab  3 \tab ("box") \tab : \tab Boxcar (Rectangular, Dirichlet)\cr
                     \tab  4 \tab ("han") \tab : \tab Hanning \cr
                     \tab  5 \tab ("ham") \tab : \tab Hamming \cr
                     \tab  6 \tab ("bla") \tab : \tab Blackman 
                     } 
                 Default: 1 = "bar".
                 }
  \item{window.size.t}{size of the window used for smoothing in time direction in units of 1/dt. Default: 5, which together with dt=1 defines a window of length 5*(1/dt) = 5. Windows of even-numbered sizes are extended by 1.}
  \item{window.size.s}{size of the window used for smoothing in scale direction in units of 1/dj. Default: 1/4, which together with dj=1/20 defines a window of length (1/4)*(1/dj) = 5. Windows of even-numbered sizes are extended by 1.}                   
  \item{make.pval}{Compute p-values? Logical. Default: TRUE.}
  \item{method}{the method of generating surrogate time series, select from:
                \tabular{rlll}{
                    \tab "white.noise"  \tab : \tab white noise \cr
                    \tab "shuffle"      \tab : \tab shuffling the given time series \cr
                    \tab "Fourier.rand" \tab : \tab time series with a similar spectrum \cr
                    \tab "AR"           \tab : \tab AR(p) \cr  
                    \tab "ARIMA"        \tab : \tab ARIMA(p,0,q)
                    }
                Default: "white.noise"
                } 
  \item{params}{a list of assignments between methods (AR, and ARIMA) and lists of parameter values applying to surrogates. Default: NULL. 
  
                Default which includes:
  
                AR: AR = list(p=1), where:
                
                \tabular{rlll}{
                    \tab p \tab : \tab AR order
                    }
                
                ARIMA: ARIMA = list(p=1, q=1, include.mean=T, sd.fac=1, trim = F, trim.prop = 0.01), where:
                
                \tabular{rlll}{  
                    \tab  p            \tab : \tab  AR order \cr
                    \tab  q            \tab : \tab  MA order \cr
                    \tab  include.mean \tab : \tab  Include a mean/intercept term? \cr
                    \tab  sd.fac       \tab : \tab  magnification factor to boost the \cr 
                    \tab               \tab   \tab  residual standard deviation \cr  
                    \tab  trim         \tab : \tab  Simulate trimmed data? \cr
                    \tab  trim.prop    \tab : \tab  high/low trimming proportion
                    }
                        
                }                         
  \item{n.sim}{number of simulations. Default: 100.}
  \item{verbose}{Print verbose output on the screen? Logical. Default: TRUE.}
}
% \details{ }
\value{A list of class \code{analyze.coherency} with the following elements:
  \item{series}{a data frame with the following columns
  
            \tabular{rllll}{
                     \tab date                 \tab : \tab the calendar date \cr
                     \tab                      \tab   \tab (as given in \code{my.data}) \cr
                     \tab <x>, <y>             \tab : \tab the two series which have been analyzed \cr
                     \tab                      \tab : \tab (detrended, if loess.span != 0; \cr
                     \tab                      \tab   \tab original names retained) \cr
                     \tab <x>.trend, <y>.trend \tab : \tab the two trend series \cr
                     \tab                      \tab   \tab (included if loess.span != 0) \cr
                     } 
                     
             Row names are resumed from \code{my.data}, and so are dates which were given as rownames.        
  }    
  \item{loess.span}{parameter alpha in loess controlling the degree of time series smoothing if the time series were detrended; no detrending if loess.span=0.}
  \item{dt}{time resolution, i.e. sampling resolution on time domain, 1/dt = number of intervals per time step.}
  \item{dj}{frequency resolution, i.e. sampling resolution on frequency domain, 1/dj = number of suboctaves (voices per octave).}
  \item{Wave.xy}{(complex-valued) cross-wavelet transform (analogous to Fourier cross-frequency spectrum, and to the covariance in statistics)} 
  \item{Angle}{phase difference, i.e. phase lead of <x> over <y> (= phase.x-phase.y)}
  \item{sWave.xy}{smoothed (complex-valued) cross-wavelet transform}
  \item{sAngle}{phase difference, i.e. phase lead of <x> over <y>, affected by smoothing} 
  \item{Power.xy}{cross-wavelet power (analogous to Fourier cross-frequency power spectrum)}
  \item{Power.xy.avg}{average cross-wavelet power in the frequency domain (averages over time)} 
  \item{Power.xy.pval}{p-values of cross-wavelet power}
  \item{Power.xy.avg.pval}{p-values of average cross-wavelet power}
  \item{Coherency}{the (complex-valued) wavelet coherency of series <x> over series <y> in the time/frequency domain, affected by smoothing (analogous to Fourier coherency, and to the coefficient of correlation in statistics)}       
  \item{Coherence}{wavelet coherence (analogous to Fourier coherence, and to the coefficient of determination in statistics (affected by smoothing)}
  \item{Coherence.avg}{average wavelet coherence in the frequency domain (averages across time)}
  \item{Coherence.pval}{p-values of wavelet coherence}
  \item{Coherence.avg.pval}{p-values of average wavelet coherence}
  \item{Wave.x, Wave.y}{(complex-valued) wavelet transforms of series <x> and <y>}
  \item{Phase.x, Phase.y}{phases of series <x> and <y>}
  \item{Ampl.x, Ampl.y}{amplitudes of series <x> and <y>}
  \item{Power.x, Power.y}{wavelet power of series <x> and <y>}
  \item{Power.x.avg, Power.y.avg}{average wavelet power of series <x> and <y>, averages across time} 
  \item{Power.x.pval, Power.y.pval}{p-values of wavelet power of series <x> and <y>}
  \item{Power.x.avg.pval, Power.y.avg.pval}{p-values of average wavelet power of series <x> and <y>}
  \item{sPower.x, sPower.y}{smoothed wavelet power of series <x> and <y>} 
  \item{Ridge.xy}{ridge of cross-wavelet power, in the form of a 0-1 matrix: columns correspond to \code{dt} steps, rows correspond to \code{dj} steps whose numerical values are given in \code{Period}}
  \item{Ridge.co}{ridge of wavelet coherence}
  \item{Ridge.x, Ridge.y}{power ridges of series <x> and <y>}
  \item{Period}{the Fourier periods (in time units)}
  \item{Scale}{the scales}
  \item{nc}{number of columns/time steps}
  \item{nr}{number of rows/scales/Fourier periods}
  \item{coi.1, coi.2}{borders of the region where the wavelet transforms are not influenced by edge effects (cone of influence)}
  \item{axis.1}{tick levels corresponding to time steps}
  \item{axis.2}{tick levels corresponding to Fourier periods (= log2(Period))}
}
\references{ 

Aguiar-Conraria L., and Soares M.J., 2011. 
Business cycle synchronization and the Euro: A wavelet analysis. 
Journal of Macroeconomics 33 (3), 477--489.

Aguiar-Conraria L., and Soares M.J., 2011. 
The Continuous Wavelet Transform: A Primer. 
NIPE Working Paper Series 16/2011.

Aguiar-Conraria L., and Soares M.J., 2012. \code{GWPackage}. Available at \url{http://sites.google.com/site/aguiarconraria/joanasoares-wavelets}; accessed September 4, 2013.

Cazelles B., Chavez M., Berteaux, D., Menard F., Vik J.O., Jenouvrier S., and Stenseth N.C., 2008.
Wavelet analysis of ecological time series.
Oecologia 156, 287--304.

Liu P.C., 1994. 
Wavelet spectrum analysis and ocean wind waves. 
In: Foufoula-Georgiou E., and Kumar P., (eds.), Wavelets in Geophysics, Academic Press, San Diego, 151--166.

Tian, H., and Cazelles, B., 2012. \code{WaveletCo}. Available at \url{http://cran.r-project.org/src/contrib/Archive/WaveletCo/}, archived April 2013; accessed July 26, 2013.

Torrence C., and Compo G.P., 1998. 
A practical guide to wavelet analysis.
Bulletin of the American Meteorological Society 79 (1), 61--78.

Veleda D., Montagne R., and Araujo M., 2012. 
Cross-Wavelet Bias Corrected by Normalizing Scales. 
Journal of Atmospheric and Oceanic Technology 29, 1401--1408.
}
\author{Angi Roesch and Harald Schmidbauer; credits are also due to Huidong Tian, Bernard Cazelles, Luis Aguiar-Conraria, and Maria Joana Soares.}
% \note{ }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ 
\code{\link{wc.avg}}, \code{\link{wc.image}}, \code{\link{wc.sel.phases}}, \code{\link{wc.phasediff.image}} 
}
\examples{
\dontrun{
## The following example is adopted from Veleda et al, 2012:

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
        ts(my.data$y, start=0, frequency=24), 
        type="l", col=1:2, 
        xlab="time (days)", ylab="hourly data",
        main="a series of hourly data with periods of 1, 2, 4, 8, 16, 32, 64, and 128 days", 
        sub="(different amplitudes at periods 16 and 64)")
legend("topright", legend=c("x","y"), col=1:2, lty=1)

## computation of cross-wavelet power and wavelet coherence:
my.wc = analyze.coherency(my.data, c("x","y"), loess.span=0, 
                          dt=1/24, dj=1/20,
                          window.size.t=1, window.size.s=1/2, 
                          lowerPeriod=1/4, 
                          make.pval=T, n.sim=10)

## plot of cross-wavelet power (with color breakpoints according to quantiles):
wc.image(my.wc, timelab="time (days)", periodlab="period (days)", 
         main="cross-wavelet power",
         legend.params=list(lab="cross-wavelet power levels"))
## plot of average cross-wavelet power:
wc.avg(my.wc)

## plot of wavelet coherence (with color breakpoints according to quantiles):
wc.image(my.wc, which.image="wc", timelab="time (days)", periodlab="period (days)", 
         main="wavelet coherence", 
         legend.params=list(lab="wavelet coherence levels", lab.line=3.5, label.digits=3))
## plot of average coherence:
wc.avg(my.wc, which.avg="wc", legend.coords="topleft")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }