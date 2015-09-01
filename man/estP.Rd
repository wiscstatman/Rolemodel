\name{estP}
\alias{estP}

\title{Estimate p using \code{mgsa}}
\description{
Uses the \code{mgsa} library to estimate \code{pi} for input to the
\code{ILP}, \code{sequentialRM}, \code{bp} functions.
}
\usage{
estP(I, y, alpha, gamma)
}

\arguments{
  \item{I}{incidence matrix, genes as the rows, gene sets as columns}
  \item{y}{binary vector, with genes as names, must match rownames(I)}
  \item{alpha}{numeric, false positive rate in role model, must be
    between 0 and 1}
  \item{gamma}{numeric, true positive rate in role model, must be
  between 0 and 1}
}

\details{
See \code{ILP}
}
\value{numeric, estimate of \code{p} for input to other Rolemodel functions}

\author{
Aimee Teo Broman, Zhishi Wang
}

\seealso{
  \code{mgsa}
  \code{ILP}
  \code{sequentialRM}
  \code{bp}
}

\examples{
data(t2d)

## estP(t2d$I,t2d$y,alpha=0.00019,gamma=0.02279)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
