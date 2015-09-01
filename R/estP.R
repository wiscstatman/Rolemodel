
estP <- function(I,y,alpha,gamma){
    Imgsa <- lapply(1:ncol(I), function(j)
                    as.character(which(I[,j]==1)))
    Imgsa <- new("MgsaSets",sets=Imgsa)
    fit <- mgsa(as.character(which(y==1)),
                Imgsa,alpha=alpha,beta=1-gamma,restarts=20)
    fit@pPost$value[which.max(fit@pPost$estimate)]
}
