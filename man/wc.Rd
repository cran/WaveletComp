\name{wc}
\alias{wc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cross-wavelet transformation, wavelet coherence computation, and a simulation algorithm}
\description{
This function provides Morlet cross-wavelet transformation results of the given two time series, 
performed within the lower-order functions \code{WaveletCoherency} and \code{WaveletTransform} 
subject to criteria concerning the time and frequency resolution, an (optional) lower and/or upper Fourier period, 
and a variety of filtering methods for the coherence computation.
It performs a simulation algorithm to test against a specified alternative hypothesis, 
which can be chosen from a variety of options, and provides p-values. 
The selected model will be fitted to the data and simulated according to estimated parameters 
in order to provide surrogate time series.
This function is called by function \code{analyze.coherency}.

The name and parts of the layout were inspired by a similar function developed by 
Huidong Tian and Bernard Cazelles (archived R package \code{WaveletCo}). 
The major part of the code for the computation of the cone of influence and the code for 
Fourier-randomized surrogate time series have been adopted from Huidong Tian.
The implementation of a choice of filtering windows for the computation of
the wavelet coherence was inspired by Luis Aguiar-Conraria and Maria Joana Soares (\code{GWPackage}).
}
\usage{
wc(x, y, start = 1, dt = 1, dj = 1/20, 
   lowerPeriod = 2*dt, upperPeriod = floor(length(x)/3)*dt, 
   window.type.t = 1, window.type.s = 1, 
   window.size.t = 5, window.size.s = 1/4,
   make.pval = TRUE, method = "white.noise", params = NULL, 
   n.sim = 100, save.sim = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the time series x to be analyzed}
  \item{y}{the time series y to be analyzed (of the same length as x}
  %%%%%%%%%%%%%%%%%
  \item{start}{starting point in time (for the computation of the cone of influence).
  
  Default: \code{start = 1}.}
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
  \item{make.pval}{Compute p-values? Logical. Default: \code{TRUE}.}
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
  \item{n.sim}{number of simulations. Default: \code{100}.}
  %%%%%%%%%%%%%%%%%
  \item{save.sim}{Shall simulations be saved on the output list? Logical.\cr Default: \code{FALSE}.}
}
% \details{ }
\value{A list with the following elements:
  %%%%%%%%%%%%%%%%%
  \item{Wave.xy}{(complex-valued) cross-wavelet transform (analogous to Fourier cross-frequency spectrum, 
  and to the covariance in statistics)} 
  \item{Angle}{phase difference, i.e. phase lead of x over y (= \code{phase.x}-\code{phase.y})}
  \item{sWave.xy}{smoothed (complex-valued) cross-wavelet transform}
  \item{sAngle}{phase difference, i.e. phase lead of x over y, affected by smoothing} 
  %%%%%%%%%%%%%%%%%
  \item{Power.xy}{cross-wavelet power (analogous to Fourier cross-frequency power spectrum)}
  \item{Power.xy.avg}{average cross-wavelet power in the frequency domain (averages over time)} 
  \item{Power.xy.pval}{p-values of cross-wavelet power}
  \item{Power.xy.avg.pval}{p-values of average cross-wavelet power}
  %%%%%%%%%%%%%%%%%
  \item{Coherency}{(complex-valued) wavelet coherency of series x over series y in the time/frequency domain,
  affected by smoothing (analogous to Fourier coherency, and to the coefficient of correlation in statistics)}       
  \item{Coherence}{wavelet coherence (analogous to Fourier coherence, and to the coefficient of determination 
  in statistics (affected by smoothing)}
  \item{Coherence.avg}{average wavelet coherence in the frequency domain (averages across time)}
  \item{Coherence.pval}{p-values of wavelet coherence}
  \item{Coherence.avg.pval}{p-values of average wavelet coherence}
  %%%%%%%%%%%%%%%%%
  \item{Wave.x, Wave.y}{(complex-valued) wavelet transforms of series x and y}
  \item{Phase.x, Phase.y}{phases of series x and y}
  \item{Ampl.x, Ampl.y}{amplitudes of series x and y}
  \item{Power.x, Power.y}{wavelet power of series x and y}
  \item{Power.x.avg, Power.y.avg}{average wavelet power of series x and y, averages across time} 
  \item{Power.x.pval, Power.y.pval}{p-values of wavelet power of series x and y}
  \item{Power.x.avg.pval, Power.y.avg.pval}{p-values of average wavelet power of series x and y}
  \item{sPower.x, sPower.y}{smoothed wavelet power of series x and y}
  %%%%%%%%%%%%%%%%%
  \item{Period}{the Fourier periods 
  (measured in time units determined by \code{dt}, see the explanations concerning \code{dt})}
  \item{Scale}{the scales (the Fourier periods divided by the Fourier factor)}
  %%%%%%%%%%%%%%%%%
  \item{coi.1, coi.2}{borders of the region where the wavelet transforms are not influenced by edge effects (cone of influence). 
  The coordinates of the borders are expressed in terms of internal axes \code{axis.1} and \code{axis.2}.} 
  %%%%%%%%%%%%%%%%%
  \item{nc}{number of columns = number of observations = number of observation epochs; "epoch" meaning point in time}
  \item{nr}{number of rows = number of scales (Fourier periods)}
  %%%%%%%%%%%%%%%%%
  \item{axis.1}{tick levels corresponding to the time steps used for (cross-)wavelet transformation: \code{1, 1+dt, 1+2dt, ...}.
  The default time axis in plot functions provided by \code{WaveletComp} is determined by observation epochs, however; "epoch" meaning point in time. }
  \item{axis.2}{tick levels corresponding to the log of Fourier periods: \code{log2(Period)}. This determines the period axis in plot functions provided by \code{WaveletComp}.}
  %%%%%%%%%%%%%%%%%
  \item{series.sim}{a data frame of the series simulated as surrogates for the (detrended) time series 
  (if both \code{make.pval = TRUE} and \code{save.sim = TRUE}.)}
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
\author{ Angi Roesch and Harald Schmidbauer; credits are also due to Huidong Tian, Bernard Cazelles,
Luis Aguiar-Conraria, and Maria Joana Soares.}
% \note{ }

%% ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{ }
% \examples{ }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }