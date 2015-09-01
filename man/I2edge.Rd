\name{I2edge}
\alias{I2edge}

\title{Create edge data frame from gene - geneset indicator matrix}
\description{Used for input to \code{bp} function, \code{edge} argument}

\usage{
I2edge(I)
}

\arguments{
  \item{I}{indicator matrix, rows are genes and columns are genesets}
}
\value{
data frame of edges:
\item{column 1}{gene sets}
\item{column 2}{genes}
}

\author{
Aimee Teo Broman
}

\seealso{
\code{gsEdge}
}
\examples{
data(t2d)

## eout <- I2edge(t2d$I)
}