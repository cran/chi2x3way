QQplot <-
function(nsample=100, yobs,nameC, taudf){
#-----------------------------------------------------------------------------------------------------------
# graphical results QQplots to check the distributin of the Marcorchino index using different theoretical probabilities 
#  observed against theoretical points
# checking the DISTRIBUTION of Marcotorchino index and its terms under Scenario 1 and 2
#----------------------------------------------------------------------------
dev.new()
x=rep(1:nsample)
x=as.matrix(t(x)/(nsample+1))
yy=as.matrix(qchisq(x, df=taudf))
zm=max(yobs)
zmm=max(yy)
zmax=max(zm,zmm)
xlimite=c(0,1.1*zmax)
ylimite=c(0,1.1*zmax)
chi_scores<-qchisq(ppoints(nsample), df = taudf)
qqplot(chi_scores, yobs,
xlab=c(nameC,taudf), ylab="Sample Quantiles ",pch=".",cex=2,col="blue",mar = c(5, 4, 4, 2) + 0.1,xlim=xlimite,ylim=ylimite)
#qqline(yobs, distribution = function(p) qchisq(p, df = taudf),col=2)
abline(a=0,b=1,col="red")
points(quantile(chi_scores,c(.01,.99)),
quantile(chi_scores,c(.01,.99)),cex=1,bg="green",pch=21)
abline(v=quantile(chi_scores,c(.01,.99)),
col="green",lty=2)
points(quantile(chi_scores,c(.05,.95)),
quantile(chi_scores,c(.05,.95)),cex=1,bg="red",pch=21)
abline(v=quantile(chi_scores,c(.05,.95)),
col="red",lty=2)
#list(chi_scores=chi_scores)
#--------------------------------------------------------------------------------------
}
