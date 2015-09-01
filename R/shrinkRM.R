shrinkRM <- function(I, y, alpha, gamma, p){

  ## Check validity of arguments
  if(p <= 0 | p >= 1)
    stop('p should be between 0 and 1')
  if(alpha <= 0 | alpha >= 1)
    stop('alpha should be between 0 and 1')
  if(gamma <= 0 | gamma >= 1)
    stop('gamma should be between 0 and 1')
  if(alpha > gamma)
    stop('alpha should be less than gamma')
  ## Check the validity of I and y
  if(!all(I %in% c(0,1)))
   stop('I should only consisit of 0 and 1')
  if(!all(y %in% c(0,1)))
   stop('y should only consisit of 0 and 1')

  
  c1 <- log(p) - log(1-p)
  c2 <- log(1-gamma) - log(1-alpha)
  c3 <- log(gamma) - log(alpha)
  Wstar <- apply(X = I, MARGIN = 2, FUN = function(x) (c1+ sum(x*y)*c3) < 0)
  
  if(sum(Wstar)==0)
  {
    newI <- I
    newy <- y
  } 
  else
  {
    if(sum(Wstar) == 1)
        indict <- any(sapply(which(I[, Wstar] == 1), function(t) sum(I[t,]*Wstar) == sum(I[t,])))
    else
        indict <- apply(X = I[, Wstar], MARGIN = 2, FUN = function(x) any(sapply(which(x == 1), function(t) sum(I[t,]*Wstar) == sum(I[t,]))))
        ##if indict[i] is TRUE, then the GO term corresponding to i must be 0 in the optimal solution of ILP. 
    Wstar[Wstar] <- indict
    W1 <- !Wstar
    if(sum(W1)==0)
    {
      print('Attention: all wholes and parts are zeros (inactive) in the optiaml solution!') ### Errors will be reported if using stop().
      return() ## NULL
    }
    if(sum(W1)==1)
    {
      newI <- as.matrix(I[,W1][I[,W1]>0])
      colnames(newI) <- colnames(I)[W1]
      newy <- y[I[,W1]>0]
    }
    if(sum(W1)>1)
    {
      P1 <- (apply(X = I[,W1], MARGIN = 1, FUN = sum) > 0)
      newI <- I[P1, W1]
      newy <- y[P1]
    }
  }
  return(list(newI = newI, newy = newy))
}
