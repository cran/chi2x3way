tau3scen1<-
function (X,pi=rep(1/dim(X)[[1]],dim(X)[[1]]),pj=rep(1/dim(X)[[2]],dim(X)[[2]]),
pk=rep(1/dim(X)[[3]],dim(X)[[3]]), digits = 3) 
{
#-------------------------------------------------------------
#prescribed probabilities-under scenario 1
#X three-way contingency table
# pi: row theoretical probabilities
# pj: column theoretical probabilities
# pk: tube theoretical probabilities
#----------------------------------------------------------------
    nn <- dim(X)
    ni <- nn[1]
    nj <- nn[2]
    nk <- nn[3]
    n <- sum(X)
    ui <- rep(1, ni)
 uj <- rep(1, nj)
   uk <- rep(1, nk)
    p3 <- X/n
    if (length(dim(X)) != 3) 
   stop("X is not a three-way table \n")
p1jk <- ui %o% pj %o% pk
fijk <- pi %o% pj %o% pk
    devt <- 1 - sum(pi^2)
#browser()
    tau3 <- sum((p3 - fijk)^2/p1jk)
itau3 <- tau3/devt
fi<- apply(p3,1,sum)
fj<- apply(p3,2,sum)
fk<- apply(p3,3,sum)
        p1j <- ui %o% pj
        p1k <- ui %o% pk
        p2ij <- pi %o% pj
        p2ik <- pi %o% pk
    khiI <- sum((fi-pi)^2) 
    khiJ <- 1/ni*sum(((fj-pj) ^2)/pj)
    khiK <-  1/ni*sum(((fk-pk) ^2)/pk)
#two-way effects
    pij <- apply(p3, c(1, 2), sum)
    pik <- apply(p3, c(1, 3), sum)
    pjk <- apply(p3, c(2, 3), sum)
    p2jk <- pj %o% pk
 #   tauij <- sum((pij-(fi %*% t(pj))-t(fj%*% t(pi)) +pi%*%t(pj))^2/p1j)
#   tauik <-  sum((pik-(fi %*% t(pk))-t(fk%*%t(pi)) + pi%*%t(pk))^2/p1k)
#browser()
   tauij <- sum((pij-(fi %*% t(pj))-t(as.matrix(fj)%*%rep(1/ni,ni)) +rep(1/ni,ni)%*%t(pj))^2/p1j)
   tauik <-  sum((pik-(fi %*% t(pk))-t(as.matrix(fk)%*%rep(1/ni,ni)) + rep(1/ni,ni)%*%t(pk))^2/p1k)
   khjk <- (1/ni)* sum((pjk-(fj %*% t(pk))-t(fk%*%t(pj)) + pj%*%t(pk))^2/p2jk)
#browser()
khin3try<-sum((p3-(pij %o% (pk))-aperm(pik%o%pj,c(1,3,2)) - aperm(pjk%o%rep(1/ni,ni),c(3,1,2))+fi%o%pj%o%pk+rep(1/ni,ni)%o%fj%o%pk+rep(1/ni,ni)%o%pj%o%fk-rep(1/ni,ni)%o%pj%o%pk  )^2/p1jk)
#cat("direct trivariate computation \n")
#print(khin3try)
#print(khin3try/devt)
#print((n-1)*(ni-1)*khin3try/devt)
#----------------------------------------------------------------------------------------------
ikhiI<-khiI/devt
ikhiJ<-khiJ/devt
ikhiK<-khiK/devt
itauij <- tauij/devt
    itauik <- tauik/devt
    ikhjk <- khjk/devt
    khin3 <- tau3 - tauij - tauik - khjk- khiI- khiJ - khiK
    ikhin3 <- khin3/devt
  #  cat("Numerator Values of partial and total indices\n")
    nom <- c("I", "J", "K","IJ", "IK", "JK", "IJK", "M")
     dres <- (ni - 1) * (nj - 1) * (nk - 1)
    di <- (ni - 1) 
    dj <- (nj - 1)
    dk <- (nk - 1)
    dij <- (ni - 1) * (nj - 1)
    dik <- (ni - 1) * (nk - 1)
    djk <- (nj - 1) * (nk - 1)
if ((pi==fi)&&(pj==fj)&&(pk==fk)) {
dtot<-dij+dik+djk+dres
di<-0
dj<-0
dk<-0
}
else{    dtot <- ni*nj*nk-1}
#   dtot <- ni*nj*nk-1
    Ci <- (n - 1) * (ni - 1) * ikhiI
  Cj <- (n - 1) * (ni - 1) * ikhiJ
  Ck <- (n - 1) * (ni - 1) * ikhiK
    Cij <- (n - 1) * (ni - 1) * itauij
    Cik <- (n - 1) * (ni - 1) * itauik
    Cjk <- (n - 1) * (ni - 1) * ikhjk
 Cijk <- (n - 1) * (ni - 1) * ikhin3
CM <- (n - 1) * (ni - 1) * itau3
    x <- c(khiI,khiJ,khiK,tauij, tauik, khjk, khin3, tau3)
    y <- (100 * x)/tau3
   zz <- c(ikhiI,ikhiJ,ikhiK,itauij, itauik, ikhjk, ikhin3, itau3)
    zz2 <- c( Ci, Cj, Ck,Cij, Cik, Cjk, Cijk, CM)
    zz3 <- c(di,dj,dk,dij, dik, djk, dres, dtot)
pvalue= 1 - pchisq(zz2, zz3)
    z <- rbind(x, zz, y, zz2, zz3,pvalue)
    nomr <- c("Numerator ", "Index", 
        "% of Inertia", "C-statistic", "df","p-value")
    dimnames(z) <- list(nomr, nom)
    z=round(z, digits = digits)
 list(z=z)
}
