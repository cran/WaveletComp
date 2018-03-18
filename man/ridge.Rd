\name{ridge}
\alias{ridge}

\title{Ridge of a wavelet power or coherence spectrum}

\description{
It finds the ridge of a (cross-)wavelet power or coherence spectrum, subject to criteria concerning:
the band of neighboring values to be included in the search area for local maximization, 
and a scale factor applied to the global maximum defining a minimum level constraint.

The basic concept of ridge determination builds on ideas developed by Huidong Tian and Bernard Cazelles 
(archived R package \code{WaveletCo}).
This code is an efficient implementation of this concept, building heavily on the use of matrices
in order to minimize computation time in R.

}
\usage{
ridge(wavelet.spectrum, band = 5, scale.factor = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wavelet.spectrum}{(cross-)wavelet power or coherence spectrum in the time/frequency domain}
  \item{band}{lower/upper band of neighboring values in the frequency domain used in the search for local maxima.\cr 
  Default: \code{5}.}
  \item{scale.factor}{a scale factor applied to the global maximum defining a minimum level constraint. \cr
  Default: \code{0.1}.}
}
% \details{ }
\value{
A matrix of \code{0}s and \code{1}s is returned which determines the course of (cross-)wavelet power or coherence ridge in the time/frequency domain.
}
\references{
Tian, H., and Cazelles, B., 2012. \code{WaveletCo}. 
Available at \url{https://cran.r-project.org/src/contrib/Archive/WaveletCo/}, archived April 2013; accessed July 26, 2013.
}
\author{Angi Roesch and Harald Schmidbauer}
\note{Credits are also due to Huidong Tian and Bernard Cazelles whose basic concept of ridge determination has been adopted.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{ }
% \examples{ }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }