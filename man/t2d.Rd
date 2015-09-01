\name{t2d}
\alias{t2d}
\docType{data}
\title{Example data set}
\description{
Type 2 diabetes data set.
}
\usage{data(t2d)}
\format{
It's a list of an incidence matrix \code{I} and the corresponding gene-level data \code{y}.
The incidence matrix is a  0-1 matrix with unique row and column names,
  where rows are genes and columns are gene-sets. Gene-level data \code{y} is a 0-1 vector with the same names as the row names of \code{I}.
}
\details{
From a large-scale genome-wide association study (GWAS) involving more
than 34,000 cases and 114,000 control subjects, 77 human genes have been
implicated as affecting T2D disease susceptibility (see reference). To
assess the functional content of this gene list, we extracted 6037 gene
ontology terms, each annotating between 5 and 50 genes. These 6037 terms
annotate a total of 10,626 genes; among the 77 T2D-associated genes, 58
are in this moderately annotated class. 
}
\references{
Zhishi W., Qiuling H., Bret L. and Michael N.: A multi-functional
analyzer uses parameter constaints to improve the efficiency of
model-based gene-set analysis (2013).  

Andrew P. M. and others: Large-scale association analysis provides
insights into the genetic architecture and pathophysiology of type 2
diabetes (2012). Nature Genetics, Volume 44-9. 
}
\examples{
data(t2d)

str(t2d$I)
str(t2d$y)
}
\keyword{datasets}
