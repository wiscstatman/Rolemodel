optimalRM <- function(xxbig, xx1, xx2)
{

      
    ind1 <- logical(ncol(xx1))

    if(ncol(xx2) == 1)
    {
         for(i in 1:ncol(xx1))
         {
              ind1[i] <- sum(xxbig[,colnames(xx1)[i]]*xxbig[,colnames(xx2)]) == sum(xx1[,i])
                       }
    }
    else
    {
         for(i in 1:ncol(xx1))
         {
              ind1[i] <- any(apply(xxbig[,colnames(xx2)], 2, function(x) sum(xxbig[,colnames(xx1)[i]]*x) == sum(xx1[,i])))
                       }
    }  
    
       
    

    return(ind1)
}

