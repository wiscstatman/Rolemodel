\name{subRM}
\alias{subRM}

\title{ Trim the incidence matrix }
\description{
Trim the incidence matrix using two parameters n.low and n.up.
}

\usage{
subRM(I, n.low, n.up)
}

\arguments{
  \item{I}{The incidence 0-1 matrix with unique row and column names, where rows are parts (genes) and columns are wholes (gene-sets).}
  \item{n.low}{The lower bound of the sums of columns in the incidence matrix.}
  \item{n.up}{The upper bound of the sums of columns in the incidence matrix.}

}

\details{ Trim the size of the incidence matrix for the computation via integer linear programming (ILP). The suggested values are n.low = 5 and n.up = 50.
 
}

\value{Return the new incidence matrix after trimming.}

\references{Zhishi W., Qiuling H., Bret L. and Michael N.: A multi-functional analyzer uses parameter constaints to improve the efficiency of model-based gene-set analysis (2013). }
%\seealso{\code{subRM}, \code{sequentialRM}}
\author{ Zhishi Wang, Michael Newton and Subhrangshu Nandi.}

\examples{
data(t2d)
## Use 5 and 10 as the trimming parameters
## newI <- subRM(t2d$I, 5, 10)
}
