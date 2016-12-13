simula <-
function(I,J,K,nran=1000,pi,pj,pk){
#The following codes generate a I x J x K contingency table under the 
#independence assumption and fixed marginal probabilities  equal to those of the observed table . 
# nran nuber of individuals -sample size-
#---------------------------------------------------
Fijk<-array(0,c(I,J,K))
#---------------------------rows
rpi<-vector("numeric",I)
#pi<-runif(I)
tot<-sum(pi)
pi<-pi/tot
rpi[1]<-pi[1]
for (i in 2:I){
rpi[i]<-rpi[i-1]+pi[i]
}
#------------------------cols
rpj<-vector("numeric",J)
#pj<-runif(J)
tot<-sum(pj)
pj<-pj/tot
rpj[1]<-pj[1]
for (j in 2:J){
rpj[j]<-rpj[j-1]+pj[j]
}
#--------------------------------tub
rpk<-vector("numeric",K)
#pk<-runif(K)
tot<-sum(pk)
pk<-pk/tot
rpk[1]<-pk[1]
for (k in 2:K){
rpk[k]<-rpk[k-1]+pk[k]
}
#-------------------------------------
for (m in 1:nran){
i<-1
rr<-runif(1)
for (ii in 1: (I-1)){
if (rr>rpi[ii])  {i=ii+1}
}
j=1
rr<-runif(1)
for (jj in 1:(J-1)){
if (rr>rpj[jj]) {j=jj+1}
}
k=1
rr<-runif(1)
for (kk in 1:(K-1)){
if (rr>rpk[kk]) {k=kk+1}
}
Fijk[i,j,k]<-Fijk[i,j,k]+1
}
list(Fijk=Fijk,pi=pi,pj=pj,pk=pk)
}
