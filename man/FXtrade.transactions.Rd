\name{FXtrade.transactions}
\alias{FXtrade.transactions}
\docType{data}
\title{Series of FX trade transactions}
\description{
USD/euro FX (foreign exchange) trade: number of worldwide transactions recorded per 5-minute intervals in July 2010. The data set contains four full weekly cycles (plus three days at the beginning of July 2010), a weekly cycle lasting from Sunday, 21:00, to Friday, 20:55. The number of transactions between Friday, 21:00, and Sunday, 20:55, is 0 or close to 0. For these intervals, variable 'active' is FALSE, otherwise TRUE.

Derived from data delivered by Morning Star.
}
\usage{data("FXtrade.transactions")}
\format{A data frame of two columns:
  \tabular{rlll}{
    \tab date         \tab : \tab date and GMT time (resolution: 5 minutes), \cr
    \tab              \tab   \tab format: "\%Y-\%m-\%d \%H:\%M:\%S" (equivalently, "\%F \%T") \cr
    \tab transactions \tab : \tab number of transactions in the 5-minute interval starting with the time indicated  \cr
    \tab active       \tab : \tab trade activity indicator
    }
}
% \details{ }
\source{
Morning Star, \url{http://www.morningstar.com/}
}
% \references{ }
\examples{
data(FXtrade.transactions)
plot(as.POSIXct(FXtrade.transactions$date, format = "\%F \%T"),
     FXtrade.transactions$transactions, type = "l",
     xlab = "day", ylab = "transactions in 5-minute intervals")
}
\keyword{datasets}
