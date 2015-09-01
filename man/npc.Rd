\name{npc}
\alias{npc}
\docType{data}
\title{Example data set}
\description{
NPC data set.
}
\usage{data(npc)}
\format{
It's a list of an incidence matrix \code{I} and the corresponding gene-level data \code{y}.
The incidence matrix is a  0-1 matrix with unique row and column names,
  where rows are genes and columns are gene-sets. Gene-level data \code{y} is a 0-1 vector with the same names as the row names of \code{I}.
}
\details{
Sengupta et al. (2006) reported a genome-wide expression study of nasopharyngeal car- cinoma (NPC), with particular attention paid to unusually extensive negative associations between host genes in the NPC cells and the expression of a key gene in the infecting Epstein- Barr virus (EBV). Genes showing extreme negative Spearman correlation were entered onto the gene list; fully 438 genes (Entrez ID) constituted the 5\% FDR-controlled list. We integrated this list with GO[5:50], 5994 GO terms that hold between 5 and 50 human genes. GO[5:50] itself annotates only about half of the human genes (10293), and only 232 of the genes showing NPC-EBV association.
}
\references{
Zhishi W., Qiuling H., Bret L. and Michael N.: A multi-functional
analyzer uses parameter constaints to improve the efficiency of
model-based gene-set analysis (2013).  
 
}
\examples{
data(npc)

str(npc$I)
str(npc$y)
}
\keyword{datasets}
