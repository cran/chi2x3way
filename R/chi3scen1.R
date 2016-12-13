chi3scen1<-
function(X,pi=rep(1/dim(X)[[1]],dim(X)[[1]]), pj=rep(1/dim(X)[[2]],dim(X)[[2]]), 
pk=rep(1/dim(X)[[3]],dim(X)[[3]]), digits=3){
#prescribed probabilities-under scenario 1
#-----------------------------------------------    
nn <- dim(X)
    ni <- nn[1]
    nj <- nn[2]
    nk <- nn[3]
    n <- sum(X)
    p3 <- X/n
    if(length(dim(X)) != 3){
        stop("X is not a three-way table\n")
    }
      pijk <- pi %o% pj %o% pk
    khi3 <- n * sum((p3-pijk)^2 /pijk)
#khi3<-n*pi*pj*pk*sum((p3-1/(pi*pj*pk))^2)
#browser()
    pij <- apply(p3, c(1, 2), sum)
   pik <- apply(p3, c(1, 3), sum)
    pjk <- apply(p3, c(2, 3), sum)
    p2ij <- pi %o% pj
        p2ik <- pi %o% pk
    p2jk <- pj %o% pk
#main effect
fi<- apply(p3,1,sum)
fj<- apply(p3,2,sum)
fk<- apply(p3,3,sum)
    khi <-n* (sum((fi-pi )^2/pi) )
    khj <-n* (sum((fj-pj) ^2/pj))
    khk <- n* (sum((fk-pk) ^2/pk))
#two-way effects
    khij <- n* sum(((pij-(fi %*% t(pj))-t(fj%*% t(pi)) +pi%*%t(pj))^2)/p2ij)
   khik <- n* sum(((pik-(fi %*% t(pk))-t(fk%*%t(pi)) + pi%*%t(pk))^2)/p2ik)
   khjk <- n* sum(((pjk-(fj %*% t(pk))-t(fk%*%t(pj)) + pj%*%t(pk))^2)/p2jk)
#----------------------------------------------------------------------------------------------
khin3 <- khi3 - khi-khj-khk-khij - khik - khjk
     nomc <- c("X2I", "X2J", "X2K","X2IJ", "X2IK", "X2JK", "X2IJK", "X2")
    x <- c(khi,khj,khk,khij, khik, khjk, khin3, khi3)
    y <- (100 * x)/khi3
    dijk <- (ni-1 ) * (nj-1 ) * (nk -1)
di<- ni-1
dj<-nj-1
dk<-nk-1
    dij <- (ni - 1) * (nj - 1)
    dik <- (ni - 1) * (nk - 1)
    djk <- (nj - 1) * (nk - 1)
    dtot<-di+dj+dk+dij+dik+djk+dijk
    df <- c(di,dj,dk,dij, dik, djk, dijk, dtot)
pvalue= 1 - pchisq(x, df)
    z <- rbind(x, y, df,pvalue)
    nomr <- c("Index ", "% of Inertia", "df", "p-value")
#browser()
    dimnames(z) <- list(nomr, nomc)
    z <- round(z, digits = digits)
    list(z = z)
}
