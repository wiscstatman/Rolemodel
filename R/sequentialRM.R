sequentialRM <- function(I, y, nupstart, by = 1, alpha, gamma, p)
{

    ## Check validity of arguments
    if(p <= 0 | p >= 1)
      stop('p should be between 0 and 1')
    if(alpha <= 0 | alpha >= 1)
      stop('alpha should be between 0 and 1')
    if(gamma <= 0 | gamma >= 1)
      stop('gamma should be between 0 and 1')
    if(alpha > gamma)
      stop('alpha should be less than gamma')
 

    new <- shrinkRM(I, y, alpha, gamma, p)
    xx <- new$newI 
    ysub <- new$newy
    ### to obtain the optimal solution on a small initial incidence matrix with nlow and nupstart.
    nlow <- min(apply(xx, 2, sum))
    inixx <- subRM(xx, nlow, nupstart)
    iniy <- ysub[rownames(inixx)]
    ininew <- shrinkRM(inixx, iniy, alpha, gamma, p)
    inires <- ILP(ininew$newI, ininew$newy, alpha, gamma, p)
    onwholes <- names(which(inires$solution[1:ncol(ininew$newI)] == 1))
    nupend <- max(apply(xx, 2, sum))

    for(i in seq(nupstart + by, nupend, by))
    {
        xxbig <- subRM(xx, n.low = nlow, n.up = i ) 
        ysubbig <- ysub[rownames(xxbig)]
        xx1 <- subRM(xx, n.low = nlow, n.up = i-by)
        xx2 <- subRM(xx, n.low = i-by+1, n.up = i ) 
        ####
        ind <- optimalRM(xxbig, xx1, xx2)
        newxx <- xxbig[, c(unique(c(colnames(xx1)[ind], onwholes)), colnames(xx2))]
        newxx <- newxx[apply(newxx, 1, sum)>0,]
        newysub <- ysubbig[rownames(newxx)]

        ####
        temp <- shrinkRM(newxx, newysub, alpha, gamma, p)
        sol <- ILP(temp$newI, temp$newy, alpha, gamma, p)
        onwholes <- names(which(sol$solution[1:ncol(temp$newI)]==1))
       
        print(i)
        
    }
   
    return(list(onwholes = onwholes, sol = sol)) 
}


