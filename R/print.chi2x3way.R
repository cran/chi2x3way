print.chi2x3way <-
function(x,digits=3, ...) {
cat("\n    Data Matrix:\n")
print(x$X)
#---------------------------------------------------------------------------
if ((x$scen==1)&&(x$indextype=="chi2") ){
cat("\n    Chi2 partition under Scenario 1\n")
}
if  ((x$scen==2)&&(x$indextype=="chi2") ) 
{
cat("\n    Chi2 partition under Scenario 2\n")
}
if  ((x$scen==1)&&(x$indextype=="tauM") ) 
{
cat("\n    Partition of the Tau index of Marcotorchino under Scenario 1\n")
}
if  ((x$scen==2)&&(x$indextype=="tauM") ) 
{
cat("\n    Partition of the Tau index of Marcotorchino under Scenario 2\n")
}
#----------------------------------------------------------------------------------------------
if ((x$simulation==TRUE)&&(x$indextype=="tauM"))
{cat("Simulation reults\n")
cat("The first 10 generated contingency tables  \n")
print(x$simulaout$XG[1:10])
if  ((x$scen==1)&&(x$indextype=="tauM") ) 
{cat("The row theoretical probabilities  of all tables \n")}
else {cat("The row theoretical probabilities  of all tables are fixed equal to the observed marginal frequencies of the original table \n")}
print(x$pi, digits=digits)
cat("The observed  row margins for  the first 10 generated tables\n")
print(x$simulaout$margI[1:10])
if  ((x$scen==1)&&(x$indextype=="tauM") ) 
{cat("The column theoretical probabilities  of all tables \n")}
else {cat("The column theoretical probabilities  of all tables are fixed equal to the observed marginal frequencies of the original table \n")}
print(x$pj, digits=digits)
cat("The observed  column  margins  for  the first 10 generated tables\n")
print(x$simulaout$margJ[1:10])
if  ((x$scen==1)&&(x$indextype=="tauM") ) 
{cat("The tube theoretical probabilities  of all tables \n")}
else {cat("The tube theoretical probabilities of all tables are fixed equal to the observed marginal frequencies of the original table \n")}
print(x$pk, digits=digits)
cat("The observed  tube margins for  the first 10 generated tables\n")
print(x$simulaout$margK[1:10])
cat("The partition of the classic formula of the  Marcotorchino index for  the first 10 generated tables\n")
print(x$simulaout$CM[1:10,])
cat("The partition of the simplified formula of the  Marcotorchino index for  the first 10 generated tables\n")
print(x$simulaout$CNewM[1:10,])
cat("The partition of the Pearson's chi-squared index for  the first ten generated tables\n")
print(x$simulaout$yphi[1:10,])
cat("The degree of freedom for each term of the three-way index partition, i.e. the theoretical mean \n")
print(x$simulaout$meanTeo,digits=digits)
cat("The mean of the observed distributions of the classical CM-statistics\n")
print(x$simulaout$meanObsC,digits=digits)
cat("The mean of the observed distributions of the simplified CM-statistics\n")
print(x$simulaout$meanObsCnew,digits=digits)
cat("The variance of the theoretical distributions\n")
print(x$simulaout$meanTeo^2,,digits=digits)
cat("The variance of the observed distributions of the classical CM-statistics\n")
print(x$simulaout$varObsC,digits=digits)
cat("The variance of the observed distributions of the simplified CM-statistics\n")
print(x$simulaout$varObsCnew,digits=digits)
}#end if simula
print(x$indexparts$z)
}
