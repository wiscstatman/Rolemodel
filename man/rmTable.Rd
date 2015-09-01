\name{rmTable}
\alias{rmTable}

\title{ One single call function to run both MAP and MCMC analysis simultaneously and to output a summary table}
\description{
Given an input of gene list, this function will do an integrative analysis including both MAP and MCMC. It will also give
a summary table about the enriched gene sets.
}

\usage{
rmTable(idlist,
                      idtype=c("SYMBOL","ENTREZ","REFSEQ",
                     "ENSEMBL","ACCNUM","UNIPROT","PMID"),
                      symbol = TRUE,
                      lib, library.loc=NULL,
                      n.upp=20, n.low=5,
                      sets=c("GO","KEGG"),
                      nupstart=10, by=1, alpha=0.01, gamma=0.02, p, nburn=1000000, 
                      ngen=10000000, sub=1000, penalty=5, initial="random",
                      annotate=TRUE)
}

\arguments{
  \item{idlist}{character vector of gene ids for
  which to test enrichment of gene set categories}
  \item{idtype}{character, type of gene ids: SYMBOL, ENTREZID, REFSEQ,
  ENSEMBL, ACCNUM, UNIPROT, PMID.  Default is SYMBOL}
  \item{symbol}{logical; if TRUE, list gene symbols in set, otherwise list gene or probe ID}
  \item{lib}{character, organism library, e.g. "org.Hs.eg"}
  \item{library.loc}{character, location of library if local}
  \item{n.upp}{numeric, upper bound for number of genes in gene sets}
  \item{n.low}{numeric, lower bound for number of genes in gene sets}
  \item{sets}{character, GO or KEGG categories}

\item{...}{for other arguments, see \code{sequentialRM} and \code{bp}}

}

\details{This function actually integrates \code{sequentialRM} and \code{bp}. That is, given an input of gene list, which can be in gene symbos, Entrez ids, Ensembles, etc., this function will run both MAP and MCMC analysis. The output of MAP is ordered according to the active probability. With the help of this function, users could have total control of the parameters for the two functions \code{sequentialRM} and \code{bp}. More imporantly, this function gives a summary data frame about the active gene sets, their GO terms, set size, number of associated genes included, and so on.
}

\value{The output is a list consisting of 4 parts: rm.table, rm.sol, bp.sol and incidence.mat.
	\item{rm.table}{A data frame of active wholes (gene sets) ordered according to the posterior active probabilities. }
	\item{rm.sol}{the output of \code{sequentialRM}.}
	\item{bp.sol}{the output of \code{bp}.}
	\item{incidence.mat}{the incidence matrix and the response vector based on the input gene list}
	}

\references{Zhishi W., Qiuling H., Bret L. and Michael N.: A multi-functional analyzer uses parameter constaints to improve the efficiency of model-based gene-set analysis (2013). }
\seealso{\code{sequentialRM}, \code{bp}}
\author{ Zhishi Wang, Michael Newton and Subhrangshu Nandi.}

\examples{
data(T2D_genelist)
idlist <- idlist[-57] ## no Entrez id for "KLHDC5"

## do not run
## res <- rmTable (idlist, lib="org.Hs.eg", n.upp=20, n.low=5, nupstart=10, by=1, alpha=0.00019, gamma=0.02279, p=0.00331)
}
