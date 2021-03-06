tau3scen2new<-
function(X, digits = 3)
{
#--------------------------------------------------------------
# Numeratore del tau tridimensionale e dei tau, chi2 marginali
# X: tabella di contingenza a tre-vie oggetto tipo array
#--------------------------------------------------------------
        nn <- dim(X)
        ni <- nn[1]
        ui <- rep(1, ni)
        nj <- nn[2]
        nk <- nn[3]
        n <- sum(X)
        p3 <- X/n
        if(length(dim(X)) != 3)
                stop("X is not a three-way table \n")
        pi <- apply(p3, 1, sum)
        pj <- apply(p3, 2, sum)
        pk <- apply(p3, 3, sum)
        pijk <- pi %o% pj %o% pk
        p1jk <- ui %o% pj %o% pk
        devt <- 1 - sum(pi^2)
        tau3 <- sum(((p3 - pijk)^2/p1jk))
        itau3 <- tau3
        pij <- apply(p3, c(1, 2), sum)
        pik <- apply(p3, c(1, 3), sum)
        pjk <- apply(p3, c(2, 3), sum)
        p1j <- ui %o% pj
        p1k <- ui %o% pk
        p2ij <- pi %o% pj
        p2ik <- pi %o% pk
        tauij <- sum(((pij - p2ij)^2/p1j))
        itauij <- tauij
        tauik <- sum(((pik - p2ik)^2/p1k))
        itauik <- tauik
        khjk <- 1/ni * (sum((pjk - (pj %o% pk))^2/(pj %o% pk)))
ikhjk<-khjk   
khin3 <- tau3 - tauij - tauik - khjk
        ikhin3 <-khin3
    #    cat("Numerator Values of partial and total indices\n")
        nom <- c("IJ", "IK", "JK", "IJK", "M")
        x <- c(tauij, tauik, khjk, khin3, tau3)
        y <- (100 * x)/tau3
        dres<-(ni-1)*(nj-1)*(nk-1)
dij<-(ni-1)*(nj-1)
dik<-(ni-1)*(nk-1)
djk<-(nj-1)*(nk-1)
dtot<-dij+dik+djk+dres
CM <- (n - 1) *(ni)* itau3  
Cij<-(n-1)*(ni)*itauij
Cik<-(n-1)*(ni)*itauik
Cjk<-(n-1)*(ni)*ikhjk
Cijk<-(n-1)*(ni)*ikhin3
zz <- c(itauij, itauik, ikhjk, ikhin3,itau3)   
zz2<-c(Cij,Cik,Cjk,Cijk,CM)
zz3<-c(dij,dik,djk,dres,dtot)
pvalue= 1 - pchisq(zz2, zz3)
     #   z1 <- rbind(zz, zz2,zz3 )
     #   dimnames(z1) <- list(c("Indices","C-Statistics","degree of freedom"),nomI)
        z <- rbind(x, zz,y,zz2,zz3,pvalue)
        nomr <- c("Numerator","Index",
"% of Inertia","CS-statistic","df","p-value")
        dimnames(z) <- list(nomr, nom)
z<-round(z,digits=digits)
        list(z=z)
}
