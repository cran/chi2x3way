tauMbootQQ <-
function(rows=3,cols=3,tubs=3,nboots=1000,nran=10000,digits=3,scen=2,
pi,pj,pk){
# REAL DATA (or one simulated array data set)
# PREAMBLE, SETTINGS
# scenario of the partition scenario 0 and 1=scen1; scenario 2=scen2 
# Number of bootstrap replicates: nboots
# set number of samples nboots=1000
# set number of individuals in each array nran=10000
# pi: row theoretical probabilities
# pj: column theoretical probabilities
# pk: tube theoretical probabilities
# checking the DISTRIBUTION of Marcotorchino index and its terms under Scenario 1 and 2
#----------------------------------------------------------------------------
tot<-rows*cols*tubs
I<-rows
J<-cols
K<-tubs
XB<-list()
XG<-list()
ristauX<-list()
ristauIX<-list()
rischiX<-list()
ristaudf<-list()
margI<-list()
margJ<-list()
margK<-list()
if (scen==1){ chi3F<-chi3scen1
tau3F<-tau3scen1
tau3FI<-tau3scen1new}
else {
chi3F<-chi3scen2
tau3F<-tau3scen2
tau3FI<-tau3scen2new}
# BOOTSTRAPPING
cont=0
g=0
for (b in 1:nboots) {
XB[[b]]<-simula(I=rows,J=cols,K=tubs,nran=nran,pi=pi,pj=pj,pk=pk)
for(i in 1: rows){
for (j in 1:cols){
for (k in 1:tubs){
if (any(XB[[b]]$pi[i]*XB[[b]]$pj[j]*XB[[b]]$pk[k]*nran<=5))  { cat("expected frequency less than 5",b,"\n")  #expected frequency n*pi*pj*pk
cont=cont+1
nboots=nboots+1
next
}
}}}
g=g+1
XG[[g]]<-XB[[b]]$Fijk # valid tables
margI[[g]]<-apply(XG[[g]]/nran,1,sum)
names(margI[[g]])<-names(XB[[g]]$pi)
margJ[[g]]<-apply(XG[[g]]/nran,2,sum)
names(margJ[[g]])<-names(XB[[g]]$pj)
margK[[g]]<-apply(XG[[g]]/nran,3,sum)
names(margK[[g]])<-names(XB[[g]]$pk)
#-------------------------------------------------------------------------------------
if (scen==1){
rischiX[[b]]<-chi3F(XB[[b]]$Fijk,pi=pi,pj=pj,pk=pk)$z[1,] # statistics chi3
ristauX[[b]]<-tau3F(XB[[b]]$Fijk,pi=pi,pj=pj,pk=pk)$z[4,] # statistics tau3
ristauIX[[b]]<-tau3FI(XB[[b]]$Fijk,pi=pi,pj=pj,pk=pk)$z[4,] # alternative formula of statistics tau3 under scen2
ristaudf[[b]]<-tau3F(XB[[b]]$Fijk,pi=pi,pj=pj,pk=pk)$z[5,] #df tau3
}
else{
rischiX[[b]]<-chi3F(XB[[b]]$Fijk)$z[1,] # statistics chi3
ristauX[[b]]<-tau3F(XB[[b]]$Fijk)$z[4,] # statistics tau3
ristauIX[[b]]<-tau3FI(XB[[b]]$Fijk)$z[4,] # alternative formula of statistics tau3 under scen2
ristaudf[[b]]<-tau3F(XB[[b]]$Fijk)$z[5,] #df tau3
#rischiX[[b]]<-chi3F(XB[[b]]$Fijk,pi=margI[[b]],pj=margJ[[b]],pk=margK[[b]])$z[1,] # statistics chi3
#ristauX[[b]]<-tau3F(XB[[b]]$Fijk,pi=margI[[b]],pj=margJ[[b]],pk=margK[[b]])$z[4,] # statistics tau3
#ristauIX[[b]]<-tau3FI(XB[[b]]$Fijk,pi=margI[[b]],pj=margJ[[b]],pk=margK[[b]])$z[4,] # alternative formula of statistics tau3 under scen2
#ristaudf[[b]]<-tau3F(XB[[b]]$Fijk,pi=margI[[b]],pj=margJ[[b]],pk=margK[[b]])$z[5,] #df tau3
}
} # end boots
#browser()
if (scen==1){ 
nameC <- c("C_I","C_J","C_K","C_IJ", "C_IK", "C_JK", "C_IJK"," C_M")
nameCnew <- c("CS_I","CS_J","CS_K","CS_IJ", "CS_IK", "CS_JK", "CS_IJK"," CS_M")
nameChi <- c("chi2_I","chi2_J","chi2_K","chi2_IJ", "chi2_IK", "chi2_JK", "chi2_IJK"," chi2")
nindex<-8
}
else{
nameC <- c("C_IJ", "C_IK", "C_JK", "C_IJK"," C_M")
nameCnew <- c("CS_IJ", "CS_IK", "CS_JK", "CS_IJK"," CS_M")
nameChi <- c("chi2_IJ", "chi2_IK", "chi2_JK", "chi2_IJK"," chi2")
nindex<-5
}
nsample=nboots-cont
ytau<-matrix(unlist(c(ristauX)),nsample,nindex,byrow=T) #Tau3-stats
ytauI<-matrix(unlist(c(ristauIX)),nsample,nindex,byrow=T) #Tau3I-stats formula alternativa sotto scen 2
ychi<-matrix(unlist(c(rischiX)),nsample,nindex,byrow=T) #Chi2-stats
taudf<-matrix(unlist(ristaudf[[1]]),1,nindex,byrow=T) #df
meanObsC<-matrix(0,1,nindex) #Tau3-stats
meanObsCnew<-matrix(0,1,nindex) #Tau3-stats
varObsC<-matrix(0,1,nindex) #Tau3-stats
varObsCnew<-matrix(0,1,nindex) #Tau3-stats
provaRSS<-matrix(0,1,nindex)
chi_scores<-matrix(0,nsample,nindex)
dimnames(taudf)<-list(NULL,nameC)
dimnames(varObsC)<-list(NULL,nameC)
dimnames(varObsCnew)<-list(NULL,nameCnew)
dimnames(meanObsC)<-list(NULL,nameC)
dimnames(meanObsCnew)<-list(NULL,nameCnew)
#browser()
for (i in 1:nindex){
chi_scores[,i]<-qchisq(ppoints(nsample), df = taudf[1,i])
## Q-Q plot for Chi^2 data against true theoretical distribution:
#---------------------------------------------------------------------------------------------
QQplot(nsample=nsample,yobs=ytau[,i],nameC=nameC[i],taudf=taudf[1,i])
#--------------------------------------------------------------------------------------
QQplot(nsample=nsample,yobs=ychi[,i],nameC=nameChi[i],taudf=taudf[1,i])
#-----------------------------------------------------------------------------------------
QQplot(nsample=nsample,yobs=ytauI[,i],nameC=nameCnew[i],taudf=taudf[1,i])
#--------------------------------------------------------------------------------------
varObsC[,i]=var(ytau[,i])
varObsCnew[,i]=var(ytauI[,i])
provaRSS[,i]=sqrt(sum((chi_scores[,i]-ytau[,i])^2))
meanObsC[,i]=mean(ytau[,i])
meanObsCnew[,i]=mean(ytauI[,i])
} #end for 
varTeo=taudf^2
ytau=round(ytau,digits=digits)
ytauI=round(ytauI,digits=digits)
ychi=round(ychi,digits=digits)
dimnames(ytau)=list(paste("sample",1:nsample,sep=""),nameC)
dimnames(ytauI)=list(paste("sample",1:nsample,sep=""),nameCnew)
dimnames(ychi)=list(paste("sample",1:nsample,sep=""),nameChi)
list(XG=XG,margI=margI,margJ=margJ,margK=margK,CM=ytau,CNewM=ytauI,yphi=ychi,meanTeo=taudf,meanObsC=meanObsC,
meanObsCnew=meanObsCnew,varTeo=varTeo,varObsC=varObsC,
varObsCnew=varObsCnew,cont=cont,provaRSS=provaRSS)
#list(XG=XG,margI=margI,margJ=margJ,margK=margK,ytau=ytau,ychi=ychi,chidf=taudf,cont=cont)
}
