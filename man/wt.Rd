\name{wt}
\alias{wt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Wavelet transformation, and a simulation algorithm}

\description{
This function provides Morlet wavelet transformation results of the given time series, 
performed within the lower-order function \code{WaveletTransform},  
subject to criteria concerning the time and frequency resolution, and an (optional) lower and/or upper Fourier period.
It performs a simulation algorithm to test against a specified alternative hypothesis, which can be chosen from a variety of options,
and provides p-values. The selected model will be fitted to the data and simulated according to estimated parameters 
in order to provide surrogate time series.
This function is called by function \code{analyze.wavelet}.

The name and parts of the layout were inspired by a similar function developed by Huidong Tian and Bernard Cazelles 
(archived R package \code{WaveletCo}).
The major part of the code for the computation of the cone of influence, and the code for Fourier-randomized surrogate time series
has been adopted from Huidong Tian.
}
\usage{
wt(x, start = 1, dt = 1, dj = 1/20, 
   lowerPeriod = 2*dt, upperPeriod = floor(length(x)*dt/3), 
   make.pval = TRUE, method = "white.noise", params = NULL, 
   n.sim = 100, save.sim = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the time series to be analyzed}
  %%%%%%%%%%%%%%%%%
  \item{start}{starting point in time (for the computation of the cone of influence).
  
  Default: \code{start = 1}.}
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
  \item{save.sim}{Should simulations be saved on the output list? Logical. Default: \code{FALSE}.}
}
% \details{ }
\value{A list with the following elements:
  %%%%%%%%%%%%%%%%%
  \item{Wave}{complex wavelet transform of the series}
  \item{Phase}{phases}
  \item{Ampl}{amplitudes}
  %%%%%%%%%%%%%%%%%
  \item{Power}{wavelet power in the time/frequency domain}
  \item{Power.avg}{vector of average wavelet power in the frequency domain (averages over time)} 
  \item{Power.pval}{p-values of wavelet power}
  \item{Power.avg.pval}{vector of p-values of average wavelet power}
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
\author{Angi Roesch and Harald Schmidbauer; credits are also due to Huidong Tian, and Bernard Cazelles}
% \note{ }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{WaveletTransform}}, \code{\link{analyze.wavelet}} } 
% \examples{ }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }