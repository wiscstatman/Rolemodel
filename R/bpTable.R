## bpout = output from bp() function

bpTable <- function(bpout,pr=0.5){
    bpord <- bpout[order(bpout$ActiveProbability,decreasing=TRUE),]
    bpord[bpord$ActiveProbability>pr,1:2]
}
