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
wc.angle(WC = WC, use.sAngle = T, p = 1, 
         which.lvl = "wp", lvl = 0, 
         which.sig = which.lvl, siglvl = 0.05, 
         col.arrow = "black")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{WC}{an object of class \code{analyze.coherency}}
  \item{use.sAngle}{Use smoothed version of phase-difference to plot arrows? Logical. Default: FALSE.}
  \item{p}{Which area should be filled with arrows displaying phase-differences?
  
           \tabular{rlll}{
                       \tab p=0 \tab : \tab area with high values of \code{which.lvl} only \cr
                       \tab     \tab   \tab (cf. \code{lvl}) \cr
                       \tab p=1 \tab : \tab area of significance of \code{which.sig} only \cr
                       \tab     \tab   \tab (cf. \code{siglvl}) \cr
                       \tab p=2 \tab : \tab area with both high values and significance \cr
                       \tab     \tab   \tab (combining p=0 and p=2)
                       }
    
           Default: 1
  }
  \item{which.lvl}{Which spectrum should be used to restrict the area of arrows according to its level? 
                           
                       \tabular{rlll}{
                                \tab "wp" \tab : \tab cross-wavelet power \cr
                                \tab "wc" \tab : \tab wavelet coherence 
                                }
                                
                        Default: "wp"

                        }
  \item{lvl}{minimum level of cross-wavelet power (or wavelet coherency) within the area of arrows (if p=0 or 2). Default: 0.}                      
  \item{which.sig}{Which spectrum and corresponding p-values should be used to restrict the area of arrows according to significance? 
                           
                       \tabular{rlll}{
                                \tab "wp" \tab : \tab cross-wavelet power (default) \cr
                                \tab "wc" \tab : \tab wavelet coherence \cr
                                \tab "wt" \tab : \tab individual wavelet power
                                }
                                
                   Default: \code{which.lvl}
  }
  \item{siglvl}{level of significance referring to \code{which.sig} (if p=1 or 2).

                      Default: 0.05
  }
  \item{col.arrow}{color of arrows. Default: "black".}
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

Tian, H., and Cazelles, B., 2012. \code{WaveletCo}. Available at \url{http://cran.r-project.org/src/contrib/Archive/WaveletCo/}, archived April 2013; accessed July 26, 2013.

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