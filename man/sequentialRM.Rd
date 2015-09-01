\name{sequentialRM}
\alias{sequentialRM}

\title{ The sequential approach }
\description{
Use the sequential approach introduced in the reference to speed up the running of integer linear programming (ILP).
}

\usage{
sequentialRM(I, y, nupstart, by = 1, alpha, gamma, p)
}

\arguments{
  \item{I}{The incidence 0-1 matrix with unique row and column names, where rows are parts (genes) and columns are wholes (gene-sets).}
  \item{y}{Gene-level 0-1 data with the same names as the row names of I.}
  \item{nupstart}{The starting upper bound used in the sequential approach.}
  \item{by}{The increment of the upper bound used in the sequential approach, default value 1.}
 \item{alpha}{The false positive rate in role model, numeric value between 0 and 1. See reference.}
  \item{gamma}{The true positive rate in role model, numeric value between 0 and 1. See reference.}
  \item{p}{The prior active probability of wholes in role model, numeric value between 0 and 1. See reference.}
}

\details{
Generally, alpha and gamma can be estimated from the gene-level data by users themselves (see reference for examples), and alpha is less than gamma. p can be estimated via R package \code{MGSA} with alpha and gamma fixed.

We first perform the ILP on an initial incidence matrix (the smaller matrix) with lower bound equal lower bound of I and upper bound nupstart; then do another ILP, making use of the results obtained from the last ILP, on the bigger incidence matrix with upper bound equal nupstart + by. This process is repeated until the original incidence matrix I is reached. The suggested value for nupstart is 10. \code{sequentialRM} is our main function to perform the ILP calculation. \code{ILP}, \code{shrinkRM} and \code{optimalRM} are all invoked in this function.
}

\value{Return a list consisting of \code{onwholes}: the active wholes, i.e., the MFA-ILP (MAP) estimate, and \code{sol}: has the same structure with the output of \code{ILP},
	\item{optimum}{the value of the objective function at the optimum}
	\item{solution}{the vector of optimal coefficients (0-1vector)}
	\item{status}{an integer with status information about the solution returned. If the control parameter canonicalize_status is set (the default) then it will return 0 for the optimal solution being found, and non-zero otherwise. If the control parameter is set to FALSE it will return the GLPK status codes.}
	}

\references{Zhishi W., Qiuling H., Bret L. and Michael N.: A multi-functional analyzer uses parameter constaints to improve the efficiency of model-based gene-set analysis (2013). }
%\seealso{}
\author{ Zhishi Wang, Michael Newton and Subhrangshu Nandi.}

\examples{
data(t2d)

## Isub <- subRM(t2d$I,5,20)
## ysub <- t2d$y[rownames(Isub)]

## set the system parameters
alpha <- 0.00019
gamma <- 0.02279
p <- 0.00331

## use the sequential approach to get the MAP estimate on a smaller
## example of type 2 diabetes 
## res <- sequentialRM(Isub, ysub, nupstart = 10, by =1, alpha, gamma, p)
}



