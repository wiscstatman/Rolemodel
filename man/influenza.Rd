\name{influenza}
\alias{influenza}
\docType{data}
\title{Example data set}
\description{
Influenza data set.
}
\usage{data(influenza)}
\format{
It's a list of an incidence matrix \code{I} and the corresponding gene-level data \code{y}.
The incidence matrix is a  0-1 matrix with unique row and column names,
  where rows are genes and columns are gene-sets. Gene-level data \code{y} is a 0-1 vector with the same names as the row names of \code{I}.
}
\details{
In a meta-analysis of four genome-wide studies of influenza virus, Hao et al. (2013) reported that 984 human genes had been detected by RNA interference as possibly being associated with viral replication. We used 3320 GO terms annotating between 5 and 50 human genes. Among the 984 influenza-involved genes, 683 are annotated to at least one of these terms.

}
\references{
Zhishi W., Qiuling H., Bret L. and Michael N.: A multi-functional
analyzer uses parameter constaints to improve the efficiency of
model-based gene-set analysis (2013).  

}
\examples{
data(influenza)

str(influenza$I)
str(influenza$y)
}
\keyword{datasets}
