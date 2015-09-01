\name{T2D_genelist}
\alias{T2D_genelist}
\docType{data}
\title{Example input gene list}
\description{
It includes 77 gene symbols which have been identified for affecting type 2 diabetes.
}
\usage{data(T2D_genelist)}
\format{
  idlist is a character vector.
}
\details{
From a large-scale genome-wide association study (GWAS) involving more
than 34,000 cases and 114,000 control subjects, 77 human genes have been
implicated as affecting T2D disease susceptibility (see reference). 
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
data(T2D_genelist)

str(idlist)
}
\keyword{datasets}
