\name{arrow}
\alias{arrow}

\title{Plot an arrow} 

\description{
It plots a polygon arrow subject to criteria concerning position, length, width, direction and color.

In plots of the cross-wavelet power (or wavelet coherency) spectrum of two time series as performed by function \code{wc.image}, 
arrows are used as a symbolic indicator of phase-differences between the series in the time/frequency domain.

This code has been developed by Huidong Tian (archived R package \code{WaveletCo}).
}
\usage{
arrow(x, y, l = 0.1, w = 0.3*l, alpha, col.arrow = 'black')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the x coordinate of the center of the arrow.}
  \item{y}{the y coordinate of the center of arrow.}
  \item{l}{the length of arrow. Default: 0.1.}
  \item{w}{the width of arrow. Default: 0.3*l.}
  \item{alpha}{the direction of the arrow, range from 0 to 2*pi.}
  \item{col.arrow}{color of arrows. Default: "black".}
}  
% \details{ }
% \value{ }
\references{
Tian, H., and Cazelles, B., 2012. \code{WaveletCo}. Available at \url{http://cran.r-project.org/src/contrib/Archive/WaveletCo/}, archived April 2013; accessed July 26, 2013.
}
\author{Code developed by Huidong Tian.}
% \note{ }

%% ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{ }
% \examples{ }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }