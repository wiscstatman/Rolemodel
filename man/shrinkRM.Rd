\name{shrinkRM}
\alias{shrinkRM}

\title{ Shrink the size of the incidence matrix }
\description{
Shrink the size of the incidence matrix by making a prejudgement about which wholes (gene-sets) and parts (genes) have to be zeros in the optimal solution of ILP. See reference for details.
}

\usage{
shrinkRM(I, y, alpha, gamma, p)
}

\arguments{
   \item{I}{The incidence 0-1 matrix with unique row and column names, where rows are parts (genes) and columns are wholes (gene-sets).}
  \item{y}{Gene-level 0-1 data with the same names as the row names of I.}
  \item{alpha}{The false positive rate in role model, numeric value between 0 and 1. See reference.}
  \item{gamma}{The true positive rate in role model, numeric value between 0 and 1. See reference.}
  \item{p}{The prior active probability of wholes in role model, numeric value between 0 and 1. See reference.}
}

\details{
Generally, alpha and gamma can be estimated from the gene-level data by users themselves (see reference for examples), and alpha is less than gamma. p can be estimated via R package \code{MGSA} with alpha and gamma fixed.

The amount of shrinkage may be dramatic, but it depends on the observed data y, the system I and system parameters alpha, gamma and p. When alpha is small and gamma is large the effects may be minimal. This function is invoked in \code{sequentialRM}.
}

\value{Return a list consisting of newI: the incidence matrix after shrinking, and newy: the corresponding part-level data.}

\references{Zhishi W., Qiuling H., Bret L. and Michael N.: A multi-functional analyzer uses parameter constaints to improve the efficiency of model-based gene-set analysis (2013). }
\seealso{\code{sequentialRM}}
\author{ Zhishi Wang, Michael Newton, Subhrangshu Nandi}

\examples{
data(t2d)
## set the system parameters
alpha <- 0.00019
gamma <- 0.02279
p <- 0.00331
\dontrun{
## shrink the matrix
## new <- shrinkRM(t2d$I, t2d$y, alpha, gamma, p)
}
}
