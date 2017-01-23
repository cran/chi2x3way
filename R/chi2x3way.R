chi2x3way<-function(X, indextype="chi2",scen = 2,
simulation=FALSE,nboots=1000,nran=1000,
pi=rep(1/dim(X)[[1]],dim(X)[[1]]),pj=rep(1/dim(X)[[2]],dim(X)[[2]]),
pk=rep(1/dim(X)[[3]],dim(X)[[3]]),digits=3){
# REAL DATA X three-way array
# chi-2 index for three-way contingency tables  and its partitions under Scenario 1 and 2
#----------------------------------------------------------------------------
if (!any(indextype==c("chi2","tauM"))) stop(paste("Indextype must be equal to chi2  or tauM"))
X<-as.array(X)
ni<-dim(X)[[1]]
nj<-dim(X)[[2]]
nk<-dim(X)[[3]]
n<-sum(X)
opi<-apply(X/n,1,sum)
opj<-apply(X/n,2,sum)
opk<-apply(X/n,3,sum)
if (scen==2){ #fun3F<-chi3new
S <- switch(indextype, "chi2" = chi3scen2(X,digits=digits),  "tauM" =tau3scen2(X,digits=digits) )
}
else {#fun3F<-chi3
S <- switch(indextype, "chi2" = chi3scen1(X,pi=pi, pj=pj, pk=pk,digits=digits),  "tauM" =tau3scen1(X,pi=pi, pj=pj, pk=pk, digits=digits) )
}
if (scen==1){ 
if (indextype=="chi2"){
nameC <- c("phi2_I","phi2_J","phi2_K","phi2_ij", "phi2_ik", "phi2_jk", "phi2_3"," phi2_tot")}
else{nameC <- c("tau_I","tau_J","tau_K","tau_ij", "tau_ik", "tau_jk", "tau_3"," tau_M")}
nindex<-8
}
else{
if (indextype=="chi2"){
nameC <- c("phi2_ij", "phi2_ik", "phi2_jk", "phi2_3"," phi2_tot")}
else{nameC <- c("tau_ij", "tau_ik", "tau_jk", "tau_3"," tau_M")}
nindex<-5
}
if ((simulation==TRUE)&&(indextype=="tauM"))
{
if (scen==2){
cat("Check the Marcotorchino index distribution using as probabilities the observed margins of your data table \n")
pi=opi
pj=opj
pk=opk
simulaout=tauMbootQQ(rows=ni,cols=nj,tubs=nk,nboots=nboots,nran=nran,digits=digits,scen=scen,pi=opi,pj=opj,pk=opk)
#Fijk=simulaout$Fijk
}
if (scen==1){
cat("The Marcotorchino index distribution using prescribed probabilities under scenario 1 \n")
simulaout=tauMbootQQ(rows=ni,cols=nj,tubs=nk,nboots=nboots,nran=nran,digits=digits,scen=scen,pi=pi,pj=pj,pk=pk)
}
}#end if simulation
if ((simulation==TRUE)&&(indextype=="chi2"))
{
cat("You do not need to check the Pearson's index distribution, the simulation study will not be performed!  \n")
simulaout<-NULL
}
if (simulation==FALSE){
simulaout<-NULL
}
cat("Results to print\n")
respart<-list(X=X,indexparts=S,indextype=indextype,simulaout=simulaout,pi=pi,pj=pj,pk=pk,scen=scen,simulation=simulation,nboots=nboots)
class(respart)<-"chi2x3way"
return(respart)
}
