#load necessary libraries
library(RcppArmadillo)
library(Rcpp)
library(MTS)
library(Matrix)
library(MASS)
library(matrixcalc)
library(statmod)

# set necessary values
data = 
N1 =       # numbers of group 1 subjects
N2 =      # numbers of group 2 subjects
P =       # numbers of ROIs
  
DataList=list() 
for(i in 1:N1){
  
}
for(i in 1:N2){
  
}
# Test, null hypothesis: the networks in the two groups are the same.
# hence the subject-specific networks are exchangeable.
# bring in Rcpp functions
sourceCpp('C:/Users/Clancy/Desktop/longitudinalFC-master/LongitInclude.cpp')
NPerm=5000
CSFit=FCanalysis(datalist=DataList, N1=N1, N2=N2, Nperms=NPerm)
CSFit$T_global
permp((1-ecdf(CSFit$T_dist_global)(CSFit$T_global))*NPerm,NPerm,N1,N2)

Q = choose(P,2) 
# multiple correlation
PLocal=rep(NA,Q)  
for(q in 1:Q){
  PLocal[q]=permp((1-ecdf(CSFit$T_dist_local[,q])(CSFit$T_local[q]))*NPerm,NPerm,N1,N2)
}
Padj=round(p.adjust(PLocal,method="BH"),4)