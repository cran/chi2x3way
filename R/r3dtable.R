r3dtable <-function(I=3,J=3,K=3,pi=rep(1/I,I),pj=rep(1/J,J),pk=rep(1/K,K),
nboots=1000,nran=10000,digits=3){
# PREAMBLE, SETTINGS
# Number of samples: nboots
# set number of samples nboots=1000
# set number of individuals in each array nran=10000
# pi: row theoretical probabilities
# pj: column theoretical probabilities
# pk: tube theoretical probabilities
#----------------------------------------------------------------------------
XB<-list()
margI<-list()
margJ<-list()
margK<-list()
cont=0
g=0
for (b in 1:nboots) {
XB[[b]]<-simula(I=I,J=J,K=K,nran=nran,pi=pi,pj=pj,pk=pk)
#for(i in 1: I){
#for (j in 1:J){
#for (k in 1:K){
#if (any(XB[[b]]$pi[i]*XB[[b]]$pj[j]*XB[[b]]$pk[k]*nran<=5))  { cat("expected frequency less than 5",b,"\n")  #expected frequency n*pi*pj*pk
#cont=cont+1
#nboots=nboots+1
#next
#}
#}}}
#g=g+1
margI[[b]]<-apply(XB[[b]]$Fijk/nran,1,sum)
names(margI[[b]])<-names(XB[[b]]$pi)
margJ[[b]]<-apply(XB[[b]]$Fijk/nran,2,sum)
names(margJ[[b]])<-names(XB[[b]]$pj)
margK[[b]]<-apply(XB[[b]]$Fijk/nran,3,sum)
names(margK[[b]])<-names(XB[[b]]$pk)
} # end boots
#-------------------------------------------------------------------------------------
list(XB=XB,margI=margI,margJ=margJ,margK=margK)
}
