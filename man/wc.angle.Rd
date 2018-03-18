\name{wc.angle}
\alias{wc.angle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot arrows into a cross-wavelet or wavelet coherency spectrum}
\description{
It adds arrows as symbolic indicators of phase-differences into 
a cross-wavelet spectrum or wavelet coherency spectrum of two time series. 

In particular, the area to be filled with
arrows can be determined in several ways: to reflect significance
(at a given level) with respect to cross-wavelet power, wavelet
coherency, or individual wavelet power, and/or to flag a
high-value region ("high" according to a given level) with respect
to cross-wavelet power (coherency values, respectively). 

There is an option to choose "smoothed" arrows (phase differences) for plotting; 
these have been computed from smoothing filters as defined in \code{analyze.coherency}.

The name and layout were inspired by a similar function developed by Huidong Tian and Bernard Cazelles 
(archived R package \code{WaveletCo}). The code for the arrow design to reflect phase differences 
has been adopted from Huidong Tian.

}
\usage{
wc.angle(WC = WC, exponent = 1, use.sAngle = TRUE, p = 1, 
         which.lvl = "wp", lvl = 0, 
         which.sig = which.lvl, siglvl = 0.05, 
         col.arrow = "black")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{WC}{an object of class \code{analyze.coherency}}
  %%%%%%%%%%%%%%%%%%%%%%%%%
  \item{exponent}{Exponent applied to cross-wavelet power (or wavelet coherence, depending on \code{which.lvl}) values before selecting the area 
  to be filled with arrows according to level, see \code{lvl} (only effective if \code{p = 0} or \code{p = 2}); the exponent should be positive.
  
  Default: \code{1}.
  }
  \item{use.sAngle}{Use smoothed version of phase difference to plot arrows? Logical.\cr Default: \code{FALSE}.}
  %%%%%%%%%%%%%%%%%%%%%%%%%
  \item{p}{Which area should be filled with arrows displaying phase differences?\cr 
  \itemize{
  \item[\code{0} :] Select the area w.r.t. power (or coherence, depending on \code{which.lvl}.\cr 
  The value of \code{lvl} is the threshold.\cr 
  Set \code{lvl = 0} to plot arrows over the entire area.
  \item[\code{1} :] Select the area w.r.t. significance as specified in \code{which.sig}.\cr 
  The value of \code{siglvl} is the threshold.
  \item[\code{2} :] Intersection of both areas above. 
  }
  Default: \code{1}.
  }  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  \item{which.lvl}{Which spectrum should be used to restrict the area of arrows according to its level? 
                           
                       \tabular{rlll}{
                                \tab \code{"wp"} \tab : \tab cross-wavelet power \cr
                                \tab \code{"wc"} \tab : \tab wavelet coherence \cr
                                }
                                
  Default: \code{"wp"}.
  }
  \item{lvl}{minimum level of cross-wavelet power (or wavelet coherence, depending on \code{which.lvl}) for arrows to be plotted 
  (only effective if \code{p = 0} or \code{p = 2}).
  
  Default: \code{0}.
  } 
  %%%%%%%%%%%%%%%%%%%%%%%%%
  \item{which.sig}{Which spectrum and corresponding p-values should be used to restrict the area of arrows according to significance? 
                       
                      \tabular{rlll}{
                                \tab \code{"wp"} \tab : \tab cross-wavelet power \cr
                                \tab \code{"wc"} \tab : \tab wavelet coherence \cr
                                \tab \code{"wt"} \tab : \tab individual wavelet power
                                }
                                
  Default: \code{which.lvl}
  }
  \item{siglvl}{level of significance referring to \code{which.sig} 
  (only effective if \code{p = 1} or \code{p = 2}).

  Default: \code{0.05}.
  }
  %%%%%%%%%%%%%%%%%%%%%%%%%
  \item{col.arrow}{arrow color. Default: "black".}
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

Tian, H., and Cazelles, B., 2012. \code{WaveletCo}. 
Available at \url{https://cran.r-project.org/src/contrib/Archive/WaveletCo/}, archived April 2013; accessed July 26, 2013.

Torrence C., and Compo G.P., 1998. 
A practical guide to wavelet analysis.
Bulletin of the American Meteorological Society 79 (1), 61--78.

Veleda D., Montagne R., and Araujo M., 2012. 
Cross-Wavelet Bias Corrected by Normalizing Scales. 
Journal of Atmospheric and Oceanic Technology 29, 1401--1408.
}
\author{Angi Roesch and Harald Schmidbauer; credits are also due to Huidong Tian, and Bernard Cazelles.}
% \note{ }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{analyze.coherency}}, \code{\link{wc.image}}, \code{\link{arrow}}
}
% \examples{ }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }