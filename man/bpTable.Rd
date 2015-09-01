\name{bpTable}
\alias{bpTable}

\title{
Wrapper for MCMC with bipartite graph
}
\description{
Ordered subset of "whole" results from \code{bp}
}
\usage{
bpTable(bpout, pr = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bpout}{data frame, output from \code{bp}}
  \item{pr}{numeric, lower bound of ActiveProbability for subsetting,
    default is pr=0.5}
}
\details{
See \code{bp}
}
\value{
data frame with "whole" results, ordered by ActiveProbability, subsetted
for ActiveProbability > pr
}
\author{
Aimee Teo Broman, Zhishi Wang
}
\seealso{
  \code{gs2edge}
  \code{bp}
}
\examples{
library(org.Hs.eg.db)
## idlist <- sample(mappedRkeys(org.Hs.egSYMBOL),1000)
## eout <- gs2edge(idlist,n.upp=30,idtype="SYMBOL",lib="org.Hs.eg")

## bpout <- bp(whole=colnames(eout$I),part=eout$y,edge=eout$edge)

## bpTable(bpout,pr=0.3)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
