subRM <- function(I, n.low, n.up)
{
    temp <- apply(I, 2, sum)
    if(sum((n.low<=temp)&(temp<=n.up)) == 1)
    {
       newI <- as.matrix(I[,(n.low<=temp)&(temp<=n.up)][I[,(n.low<=temp)&(temp<=n.up)]==1])
       colnames(newI) <- colnames(I)[(n.low<=temp)&(temp<=n.up)]
    }
    else
    {
       newI <- I[,(n.low<=temp)&(temp<=n.up)]
       newI <- newI[apply(newI,1,sum)>0,]
    }

    return(newI)
}
