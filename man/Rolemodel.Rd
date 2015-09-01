\name{Rolemodel}
\alias{Rolemodel}

\title{ One single call function to run both MAP and MCMC analysis simultaneously }
\description{
Given an input of gene list, this function will do an integrative analysis including both MAP and MCMC.
}

\usage{
Rolemodel(idlist, idtype=c("SYMBOL","ENTREZ","REFSEQ", "ENSEMBL","ACCNUM","UNIPROT","PMID"), 
lib, library.loc=NULL, n.upp=20, n.low=5, sets=c("GO","KEGG"), nupstart=10, by=1, alpha=0.01, gamma=0.02, p, nburn=1000000, 
ngen=10000000, sub=1000, penalty=5, initial="random", annotate=TRUE)
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

\item{...}{for other arguments, see \code{sequentialRM} and \code{bp}}

}

\details{This function actually integrates \code{sequentialRM} and \code{bp}. That is, given an input of gene list, which can be in gene symbos, Entrez ids, Ensembles, etc., this function will run both MAP and MCMC analysis. The output of MAP is ordered according to the active probability. With the help of this function, users could have total control of the parameters for the two functions \code{sequentialRM} and \code{bp}. 
}

\value{The output is a list consisting of three parts: setprobs, rmsol and bpsol.
	\item{setprobs}{A data frame of active wholes ordered according to the posterior active probabilities. }
	\item{rmsol}{the output of \code{sequentialRM}.}
	\item{bpsol}{the output of \code{bp}.}
	}

\references{Zhishi W., Qiuling H., Bret L. and Michael N.: A multi-functional analyzer uses parameter constaints to improve the efficiency of model-based gene-set analysis (2013). }
\seealso{\code{sequentialRM}, \code{bp}}
\author{ Zhishi Wang, Michael Newton and Subhrangshu Nandi.}

\examples{
data(T2D_genelist)
idlist <- idlist[-57] ## no Entrez id for "KLHDC5"

## res <- Rolemodel (idlist, lib="org.Hs.eg", n.upp=20, n.low=5, nupstart=10, alpha=0.00019, gamma=0.02279, p=0.00331, ## nburn=1000000, ngen=10000000, sub=1000, penalty=5, initial="random", annotate=TRUE)
}
