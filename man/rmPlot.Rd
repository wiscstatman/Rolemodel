\name{rmPlot}
\alias{rmPlot}

\title{ Display an image of enriched gene sets by MFA analysis }
\description{
Plot a waterfall plot about the gene sets according to their overlaps with the input gene list.
}

\usage{
rmPlot(I, y, onwholes, sets=c("GO", "KEGG"))
}

\arguments{
  \item{I}{The incidence 0-1 matrix with unique row and column names, where rows are parts (genes) and columns are wholes (gene-sets).}
  \item{y}{Gene-level 0-1 data with the same names as the row names of I.}
  \item{onwholes}{Character, the list of enriched gene sets}
  \item{sets}{Character, GO or KEGG categories.}
  
}

\details{
The enriched gene set having the highest overlap with the input gene list is plotted at the top row; then this gene set and genes in this set are removed from the system. The 2nd row will plot the gene set having the highest overlap with the remainder, and so on.
}


\references{Zhishi W., Qiuling H., Bret L. and Michael N.: A multi-functional analyzer uses parameter constaints to improve the efficiency of model-based gene-set analysis (2013). }
\seealso{\code{rmTable}}
\author{ Zhishi Wang, Michael Newton and Subhrangshu Nandi.}

\examples{
## do not run
## res <- rmTable (idlist, lib="org.Hs.eg", n.upp=20, n.low=5, nupstart=10, by=1, alpha=0.00019, gamma=0.02279, p=0.00331)
## rmPlot(res[[4]]$I, res[[4]]$y, res[[2]]$onwholes)
}
