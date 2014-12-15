\name{ridge}
\alias{ridge}

\title{Ridge of a wavelet power or coherence spectrum}

\description{
It finds the ridge of a wavelet power or coherence spectrum, subject to criteria concerning:
the band of neighboring values to be included in the search area for local maximization, 
and a scale factor applied to the global maximum defining a minimum level constraint.

This code is an efficient implementation of the basic concept of ridge determination by Huidong Tian and Bernard Cazelles 
(archived R package \code{WaveletCo}).
}
\usage{
ridge(wavelet.spectrum, band = 5, scale.factor = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wavelet.spectrum}{wavelet power or coherence spectrum in the time/frequency domain}
  \item{band}{lower/upper band of neighboring values in the frequency domain used in the search for local maxima.\cr Default: 5.}
  \item{scale.factor}{a scale factor applied to the global maximum defining a minimum level constraint. \cr
             Default: 0.1.}
}
% \details{ }
\value{
A matrix of zeros and ones is returned which determines the course of wavelet power or coherence ridge in the time/frequency domain.
}
\references{
Tian, H., and Cazelles, B., 2012. \code{WaveletCo}. Available at \url{http://cran.r-project.org/src/contrib/Archive/WaveletCo/}, archived April 2013; accessed July 26, 2013.
}
\author{Angi Roesch and Harald Schmidbauer}
\note{Credits are also due to Huidong Tian and Bernard Cazelles whose basic concept of ridge determination has been adopted.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{ }
% \examples{ }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }