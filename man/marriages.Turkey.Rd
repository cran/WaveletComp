\name{marriages.Turkey}
\alias{marriages.Turkey}
\docType{data}
\title{Series of monthly marriages in Turkey}
\description{
Series of monthly marriages in Turkey (1988-2013), as reported by DIE (Devlet Istatistik Enstitusu) / TUIK (Turkiye Istatistik Kurumu).
}
\usage{data("marriages.Turkey")}
\format{ A data frame of three columns:
        \tabular{rlll}{
                 \tab date          \tab : \tab end-of-month date \cr   
                 \tab               \tab   \tab format: "\%Y-\%m-\%d" \cr
                 \tab n.Sun         \tab : \tab number of Sundays in this month \cr
                 \tab marriages     \tab : \tab number of marriages in this month
                 }
}                 
% \details{ }
\source{ DIE (Devlet Istatistik Enstitusu) / TUIK (Turkiye Istatistik Kurumu) 

Jan 1988 to Dec 2000:\cr "Evlenme istatistikleri", DIE (Devlet Istatistik Enstitusu, Ankara), ISSN: 1300-1086; several issues.

Jan 2001 to Dec 2013:\cr \url{http://www.tuik.gov.tr/VeriTabanlari.do?vt_id=21&ust_id=109}; accessed Oct 1, 2014.
}
% \references{ }
\examples{
data(marriages.Turkey)
plot(ts(marriages.Turkey$marriages, frequency=12, start=c(1988,1)), type="l", 
     xlab="", ylab="", 
     main="monthly marriages in Turkey (1988-2013)")     
}
\keyword{datasets}
