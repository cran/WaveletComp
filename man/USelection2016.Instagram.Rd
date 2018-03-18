\name{USelection2016.Instagram}
\alias{USelection2016.Instagram}
\docType{data}
\title{Hourly time series of the number of candidate-related media posted on Instagram during the week before the 2016 US presidential election}

\description{One week (Sunday, 2016-10-30 23:00:00 EDT through Sunday, 2016-11-06 23:00:00 EST; 170 hours) of hourly readings of the number of media posted on Instagram, supposedly (according to hashtags with which they were annotated) positive/neutral or in opposition towards candidates Trump and Clinton. The 2016 US presidential election took place on Tuesday, 2016-11-08.
}

\usage{data("USelection2016.Instagram")}
\format{A data frame of four columns:
  \tabular{rlll}{
    \tab date          \tab : \tab hour of measurement (in EST5EDT, Eastern Time Zone), \cr
    \tab               \tab   \tab format: \code{"\%F \%T"} \cr
    \tab trump.pos     \tab : \tab number of media uploads with hashtag suggesting \cr
    \tab               \tab   \tab positive/neutral annotation to Trump \cr
    \tab clinton.pos   \tab : \tab number of media uploads with hashtag suggesting \cr
    \tab               \tab   \tab positive/neutral annotation to Clinton \cr
    \tab trump.neg     \tab : \tab number of media uploads with hashtag suggesting \cr
    \tab               \tab   \tab negative annotation to Trump \cr
    \tab clinton.neg   \tab : \tab number of media uploads with hashtag suggesting \cr
    \tab               \tab   \tab negative annotation to Clinton
    }
}
\details{
Media posted on Instagram are usually annotated with hashtags. A hashtag can be used to determine whether a candidate-related posting is positive/neutral towards a candidate (for example, #makeamericagreatagain for Trump, #hillary2016 for Clinton) or opposing a candidate (for example, #dumptrump for Trump, #neverhillary for Clinton). In this way, four hourly time series are obtained: Trump vs. Clinton, supporters vs. opponents. For further details, see the reference below. The date column also contains a label EDT (Eastern Daylight Time) or EST (Eastern Standard Time); daylight saving time ended 2016-11-06 at 1:00 a.m. when clocks were moved back to 1:00 a.m. EST. The time stamp "2016-11-06 02:00:00" therefore occurs twice, once with EDT and once with EST.
}
\source{
Hourly readings of the number of media posted on Instagram were collected using Instagram's built-in API.
}
\references{
Schmidbauer H., Roesch A., Stieler F., 2018.
The 2016 US presidential election and media on Instagram: Who was in the lead?
Computers in Human Behavior 81, 148--160. 
\doi{10.1016/j.chb.2017.11.021}  
}
\examples{

data(USelection2016.Instagram)
attach(USelection2016.Instagram)

my.date <- as.POSIXct(date, format = "\%F \%T", tz = "EST5EDT")

plot(my.date, trump.pos, type = "l", col = 1, lwd = 2,
     ylab = "number of media posted on Instagram", ylim = c(0,6e+6),
     xlab = "the week before the Election Day (Tuesday, 2016-11-08)")
lines(my.date, clinton.pos, col = 2, lwd = 2) 
lines(my.date, trump.neg, col = 3, lwd = 2)
lines(my.date, clinton.neg, col = 4, lwd = 2) 
legend("topleft", legend=names(USelection2016.Instagram[-1]), 
       lty = 1, lwd = 2, col = 1:4) 

detach(USelection2016.Instagram)

}
\keyword{datasets}
