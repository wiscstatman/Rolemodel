# Rolemodel
R package for gene-set analysis, deploying the multi-functional analyzer (MFA)

The role model is a probability model used in the context of gene set analysis to describe genome-wide binary, gene-level data, such as the vector of indicators recording whether or not each gene exhibits interesting features in some experiment.  The user supplies such a vector (equivalently, a list of genes exhibiting interesting properties, such as differential expression), and the MFA produces a list of gene sets from a knowledge base (e.g. Gene Ontology) that aims to describe the gene list. 

Read about the role model and MFA [here] 
(http://projecteuclid.org/euclid.aoas/1430226091) or [here]
(http://arxiv.org/abs/1310.6322)

MFA is an example of [multi-set analysis] 
 (http://www.annualreviews.org/eprint/A7Se8wheXD5rTtmituTk/full/10.1146/annurev-statistics-010814-020335).  
Gene sets are output not only owing to being enriched for interesting genes, but because the output taken together concisely describes the input gene list.  Compare to simpler gene set enrichment schemes, most of which calculate statistics one set at a time.

The role model is a generative probability model for gene-level data.  It associates to each set a binary *activity* variable, that is either null (0) or non-null (1).  Genes are assigned to multiple sets, owing in part to the multiple biological roles that they play in the cell.  A gene is said to be active if it is annotated to any activated set.  The MFA inference tool works out properties of the posterior distribution of activity states of sets given observed binary data on genes.

 The MFA output is effective when most of the interesting genes are in one of the output sets, when most of the uninteresting genes are not in any of the output sets, and when the output sets have relatively little overlap.   Formally, the MFA deploys a Bayesian inference strategy, reporting the maximum a posteriori (MAP) estimate of non-null gene sets as well as marginal posterior probabilities that each set is non-null.
