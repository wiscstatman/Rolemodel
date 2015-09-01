\name{gs2edge}
\alias{gs2edge}

\title{
Create Gene set matrices
}
\description{
Given a vector of genes and a gene set library, create input objects for
enrichment functions
}
\usage{
gs2edge(idlist, idtype = c("SYMBOL", "ENTREZID", "REFSEQ", "ENSEMBL", "ACCNUM", "UNIPROT", "PMID"), lib, library.loc = NULL, n.upp = 50, n.low = 5, sets = c("GO", "KEGG"))
}

\arguments{
  \item{idlist}{character vector of gene ids for
  which to test enrichment of gene set categories}
  \item{idtype}{character, type of gene ids: SYMBOL, ENTREZID, REFSEQ,
  ENSEMBL, ACCNUM, UNIPROT, PMID.  Default is SYMBOL}
  \item{lib}{character, organism library, e.g. "org.Hs.eg"}
  \item{library.loc}{character, location of library if local}
  \item{n.upp}{numeric, upper bound for number of genes in gene sets}
  \item{n.low}{numeric, lower bound for number of genes in gene sets}
  \item{sets}{character, GO or KEGG categories}
}

\value{a list:
\item{I}{Indicator matrix of genes (rows) and gene sets (columns), for
  gene sets of size n, where n.low <= n <= n.upp}
\item{y}{Indicator vector of genes in gene sets that are in 'idlist'}
\item{edge}{data.frame of edges between genes and gene sets}
}

\author{Aimee Teo Broman}

\seealso{
  \code{sequentialRM}
  \code{ILP}
  \code{bpList}
}
\examples{
library(org.Hs.eg.db)
## idlist <- sample(mappedRkeys(org.Hs.egSYMBOL),1000)

## out <- gs2edge(idlist,n.upp=30,idtype="SYMBOL",lib="org.Hs.eg")

## For input to sequentialRM:
## I=out$I
## y=out$y

## For input to bp:
## whole=colnames(out$I)
## part=out$y
## edge=out$edge

}
