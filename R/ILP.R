ILP <- function(I, y, alpha, gamma, p)
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
    ## Check the validity of I and y
   if(!all(I %in% c(0,1)))
      stop('I should only consist of 0 and 1')
   if(!all(y %in% c(0,1)))
      stop('y should only consist of 0 and 1')   
   
   m <- nrow(I)
   n <- ncol(I)
   ##prepare the matrix to identify the inequalities 
   mat <- matrix(0, m + n + sum(I), m + n)

   ##generate A_p - Z_w >= 0 for all I[p,w]=1
   count <- 0
   for (i in 1:n)
   {
      ind <- which(I[,i] == 1)
      for(j in 1:length(ind))
      {
	 mat[count+j,i] <- -1
	 mat[count+j,n+ind[j]] <- 1
      }
      count <- count + length(ind)
   }
   ##generate \sum_w Z_w - A_p >= 0 for all I[p,w]=1
   for(i in 1:m)
   {
      mat[count+i, which(I[i,] == 1)] <- 1
      mat[count+i, n+i] <- -1
   }
   count <- count + m
   ## generate \sum_p (Z_w-2A_p+2)>=1 for all I[p,w]=1
   sub_rhs <- numeric(n)
   for(i in 1:n)
   {
      mat[count+i, i] <- sum(I[,i])
      mat[count+i, n+which(I[,i] == 1)] <- -2
      sub_rhs[i] <- 1 - 2*sum(I[,i])
   }
   ##ILP solver
   dir <- rep(">=", m+n+sum(I))
   rhs <- c(rep(0,m+sum(I)), sub_rhs)
   types <- rep("B", m+n)
   max <- TRUE

   ##the objective function
   sub_obj <- (log(gamma)-log(alpha)-log(1-gamma)+log(1-alpha))*y + log(1-gamma)-log(1-alpha)
   obj <- c(rep(log(p)-log(1-p), n), sub_obj)
   
   sol <- Rglpk_solve_LP(obj, mat, dir, rhs, types = types, max = max)
   names(sol$solution) <- c(colnames(I), rownames(I)) 
   return(sol)
}
