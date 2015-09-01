# MCMC with bipartite graph
#
#  whole Vector of character strings with names of whole nodes
#  part Vector of 0's and 1's indicating active part nodes; names are the names of the part nodes
#  edge Matrix with two columns; each row indicates an edge between a whole node and a part node
#  alpha Parameter alpha (0 < alpha < gamma < 1)
#  gamma Parameter alpha (0 < alpha < gamma < 1)
#  p Parameter p (0 < p < 1)
#  nburn Number of burn-in generations
#  ngen Number of sample generations
#  sub Subsample rate for burn-in and sample files
#  penalty Penalty per illegal node to loglikelihood
#  initial Initial state (see Details)
#
# return value: data frame with "whole" results; burn+sample detail are included as an attribute,
# "samples", which is a matrix
#
# The initial argument can take one of three values:
# "inactive" - all whole nodes inactive; "random" - all
# whole nodes active with probability p, no illegal nodes; or
# "high" - all nodes with proportion of connected part nodes
# with response equal to 1 above 0.4 are active, no illegal nodes.
#
# Example:
# data(t2d)
# bp.out <- bp(whole=t2d$whole, part=t2d$part, edge=t2d$edge,
#              nburn=1000, ngen=1000, sub=100)
bp <-
function(whole, part, edge,
         alpha=0.05, gamma=0.2, p=0.01,
         nburn=10000, ngen=100000, sub=1000, penalty=2,
         initial=c("inactive", "random", "high"))
{
  initial <- match.arg(initial)

  stopifnot(alpha > 0, alpha < gamma, gamma < 1)
  stopifnot(p > 0, p < 1)
  stopifnot(nburn >=0)
  stopifnot(ngen >= 0)
  stopifnot(sub >= 0)
  stopifnot(penalty >= 0)

  stopifnot(length(whole) == length(unique(whole)))
  stopifnot(length(names(part)) == length(unique(names(part))))
  stopifnot(all(!is.na(part) & part==0 | part==1))
  stopifnot(all(edge[,1] %in% whole), all(edge[,2] %in% names(part)))

  maxnchar <- max(nchar(whole))
  blank <- paste(rep(" ", maxnchar), collapse="")

  nsaved <- ceiling(nburn/sub) + ceiling(ngen/sub)
  nsavedDbl <- 3
  nsavedInt <- 7

  z <- .C("R_bp",
          as.integer(length(whole)),
          as.character(whole),
          as.integer(length(part)),
          as.character(names(part)),
          as.integer(part),
          as.integer(nrow(edge)),
          as.character(edge[,1]),
          as.character(edge[,2]),
          as.double(alpha),
          as.double(gamma),
          as.double(p),
          as.integer(nburn),
          as.integer(ngen),
          as.integer(sub),
          as.double(penalty),
          as.character(initial),
          resultName=as.character(rep(blank, length(whole))),
          as.integer(maxnchar),
          resultProb=as.double(rep(0, length(whole))),
          resultCount=as.integer(rep(0, length(whole))),
          resultSample=as.integer(rep(0, length(whole))),
          resultDegree=as.integer(rep(0, length(whole))),
          resultResponse=as.integer(rep(0, length(whole))),
          as.integer(nsaved),
          savedDouble=as.double(rep(0, nsaved*nsavedDbl)),
          savedInt=as.integer(rep(0, nsaved*nsavedInt)),
          PACKAGE="Rolemodel")

  result <- data.frame(Name=z$resultName,
                       ActiveProbability=z$resultProb,
                       Count=z$resultCount,
                       Sample=z$resultSample,
                       Degree=z$resultDegree,
                       Response=z$resultResponse, stringsAsFactors=FALSE)

  saved <- cbind(matrix(z$savedDouble, ncol=nsavedDbl), matrix(z$savedInt, ncol=nsavedInt))
  colnames(saved) <- c("AcceptanceProb", "LogLikelihood", "LogPrior",
                       "ActiveOne", "ActiveZero", "InactiveOne", "InactiveZero", "Illegal", "ActiveWhole", "ActivePart")
  attr(result, "samples") <- saved

  result
}
