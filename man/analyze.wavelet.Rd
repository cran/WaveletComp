\name{analyze.wavelet}
\alias{analyze.wavelet}

\title{
Computation of the wavelet power spectrum of a single time series
}

\description{
The time series is selected from an input data frame by specifying
either its name or its column number. Optionally, the time series
is detrended, using \code{loess} with parameter \code{loess.span}.
Internally, the series will be further standardized before it
undergoes wavelet transformation.

The wavelet power spectrum is computed by applying the Morlet
wavelet. P-values to test the null hypothesis that a period (within lowerPeriod and upperPeriod) is irrelevant at a certain time
are calculated if desired; this is accomplished with the help of a
simulation algorithm. There is a selection of models from which to
choose the alternative hypothesis. The selected model will be fitted to
the data and simulated according to estimated parameters in order
to provide surrogate time series.

Wavelet transformation, as well as p-value computations, are
carried out by calling subroutine \code{wt}.

The name and parts of the layout of subroutine \code{wt} were inspired by a similar function developed by 
Huidong Tian and Bernard Cazelles (archived R package \code{WaveletCo}). The basic concept of the simulation algorithm, and of ridge determination build on ideas developed by these authors. The major part of the code for the computation of the cone of influence, and the code for Fourier-randomized surrogate time series
has been adopted from Huidong Tian. 

Wavelet computation, the simulation algorithm, and ridge determination build heavily on the use of matrices in order to minimize computation time in R. 

This function provides a broad variety of final as well as intermediate results 
which can be further analyzed in detail. 
}

\usage{
analyze.wavelet(my.data, my.series = 1, loess.span = 0.75, dt = 1, dj = 1/20, 
                lowerPeriod = 2*dt, upperPeriod = floor(nrow(my.data)/3)*dt, 
                make.pval = T, method = "white.noise", params = NULL,
                n.sim = 100, verbose = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{my.data}{data frame of time series (including header, and dates as rownames or as separate column named "date" if available)}
  \item{my.series}{name or column index indicating the series to be analyzed, e.g. 1, 2, "dji", "ftse". Default: 1}
  \item{loess.span}{parameter alpha in loess controlling the degree of time series smoothing, if the time series is to be be detrended; no detrending if loess.span=0. Default: 0.75}
  \item{dt}{time resolution, i.e. sampling resolution on time domain, 1/dt = number of intervals per time unit. Default: 1.}
  \item{dj}{frequency resolution, i.e. sampling resolution on frequency domain, 1/dj = number of suboctaves (voices per octave). Default: 1/20.}
  \item{lowerPeriod}{lower Fourier period (in time units) for wavelet decomposition. Default: 2*dt.}
  \item{upperPeriod}{upper Fourier period (in time units) for wavelet decomposition. Default: (floor of one third of time series length)*dt}
  \item{make.pval}{Compute p-values? Logical. Default: TRUE}
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
  \item{params}{a list of assignments between methods (AR, and ARIMA) and lists of parameter values to apply to surrogates. Default: NULL. 
  
                Default includes:
                
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
  \item{n.sim}{number of simulations. Default: 100}
  \item{verbose}{Print verbose output on the screen? Logical. Default: TRUE}
}
\details{
Wavelet transformation, as well as p-value computations, are
carried out by calling the internal function \code{wt}.
}
\value{A list of class \code{analyze.wavelet} with the following elements:
  \item{series}{a data frame with the following columns
  
            \tabular{rllll}{
                     \tab date             \tab : \tab the calendar date \cr
                     \tab                  \tab   \tab (as given in \code{my.data}) \cr
                     \tab <x>              \tab : \tab the series which has been analyzed \cr
                     \tab                  \tab : \tab (detrended, if loess.span != 0; \cr
                     \tab                  \tab   \tab  original names retained) \cr
                     \tab <x>.trend        \tab : \tab the trend series (if loess.span != 0)
                     } 
                     
             Row names are taken over from my.data, and so are dates if given as row names.   
  }
  \item{loess.span}{parameter alpha in loess which controlled the degree of time series smoothing, if the time series was detrended; detrending was omitted if loess.span=0.}
  \item{dt}{time resolution, i.e. sampling resolution on time domain, 1/dt = number of intervals per time step.}
  \item{dj}{frequency resolution, i.e. sampling resolution on frequency domain, 1/dj = number of suboctaves (voices per octave).}
  \item{Wave}{complex wavelet transform of the series}
  \item{Phase}{phases}
  \item{Ampl}{amplitudes}
  \item{Power}{wavelet power in the time/frequency domain}
  \item{Power.avg}{average wavelet power in the frequency domain (averages over time)} 
  \item{Power.pval}{p-values of wavelet power}
  \item{Power.avg.pval}{p-values of average wavelet power}
  \item{Ridge}{power ridge, in the form of a 0-1 matrix: columns correspond to \code{dt} steps, rows correspond to \code{dj} steps whose numerical values are given in \code{Period}} 
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

Tian, H., and Cazelles, B., 2012. \code{WaveletCo}. Available at \url{http://cran.r-project.org/src/contrib/Archive/WaveletCo/}, archived April 2013; accessed July 26, 2013.

Torrence C., and Compo G.P., 1998. 
A practical guide to wavelet analysis. 
Bulletin of the American Meteorological Society 79 (1), 61--78.
}
\author{Angi Roesch and Harald Schmidbauer; credits are also due to Huidong Tian, and Bernard Cazelles.}
% \note{ }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{wt.avg}}, \code{\link{wt.image}}, \code{\link{wt.sel.phases}}, \code{\link{wt.phase.image}}, \code{\link{reconstruct}} }
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

my.wt = analyze.wavelet(my.data, "x", 
                        loess.span=0, 
                        dt=1/24, dj=1/20, 
                        lowerPeriod=1/4, 
                        make.pval=T, n.sim=10)

## Plot of wavelet power spectrum (with equidistant color breakpoints):  
wt.image(my.wt, color.key="interval", legend.params=list(lab="wavelet power levels"))
## Plot of average wavelet power:
wt.avg(my.wt, siglvl=0.05, sigcol="red")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }