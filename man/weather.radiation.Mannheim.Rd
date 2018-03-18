\name{weather.radiation.Mannheim}
\alias{weather.radiation.Mannheim}
\docType{data}
\title{Series of weather data and ambient gamma dose rate readings}
\description{
Ten years (from January 2005 through December 2014; 3652 values) of daily mean temperature, relative humidity, and ambient gamma dose rate (German abbreviation: ODL) readings from Mannheim-Rheinau (Germany).

ODL data retrieved in 2015 from ODL-INFO \url{http://odlinfo.bfs.de/download.php} of the German Federal Office for Radiation Protection (in German: Bundesamt fuer Strahlenschutz, BfS), weather data from the German Meteorological Office (in German: Deutscher Wetterdienst, DWD) \url{ftp://ftp-cdc.dwd.de/pub/CDC/observations_germany/climate/daily/kl/historical/}.\cr --- We owe this example to our former student Nadiya Appelhans.
}
\usage{data("weather.radiation.Mannheim")}
\format{A data frame of four columns:
  \tabular{rlll}{
    \tab date          \tab : \tab day of measurement (in UTC, Coordinated Universal Time), \cr
    \tab               \tab   \tab format: "\%Y-\%m-\%d" \cr
    \tab temperature   \tab : \tab mean daily temperature, in degrees Celsius \cr
    \tab humidity      \tab : \tab mean daily relative humidity, in percent \cr
    \tab radiation     \tab : \tab mean daily ambient gamma dose rate (ODL), in microsieverts per hour
    }
}
\details{The ambient gamma dose rate is an equivalent dose representing the stochastic health effects of low levels of ionizing radiation on the human body. According to BfS, radioactivity is to be found everywhere in the environment. It may be of natural or artificial origin. --- We drew Mannheim randomly from a set of places for which weather and radiation data were available; there is no conspicuous radiation in Mannheim.}
\source{
The Federal Office for Radiation Protection (BfS) Germany, \url{http://www.bfs.de/EN/home/home_node.html}, and in particular \url{http://www.bfs.de/EN/topics/ion/environment/environment_node.html}

ODL-INFO \url{http://odlinfo.bfs.de/DE/service/downloadbereich.html}

DWD (Deutscher Wetterdienst) \url{https://www.dwd.de/EN/Home/home_node.html}
}
% \references{ }
\examples{
data(weather.radiation.Mannheim)

plot(as.Date(weather.radiation.Mannheim$date, tz = "UTC"),
     weather.radiation.Mannheim$radiation, type = "l",
     xlab = "day", 
     ylab = "mean daily ambient gamma dose rate (ODL) in microsieverts per hour")
}
\keyword{datasets}
