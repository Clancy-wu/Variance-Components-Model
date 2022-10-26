library(readr)
patient_files = list.files(path="./ROIsignals/patient/", pattern=".txt", full.names = TRUE)
health_files = list.files(path="./ROIsignals/health/", pattern=".txt", full.names = TRUE)

interested_roi = c(1,2,3) # defined by yourself

DataList = list()

ROI2List <- function(file_name){
    data = read_table(file=file_name, col_names=FALSE)
    data.interest = data[ ,interested_roi]
    data.interest = as.matrix(data.interest)
    return(data.interest)
}
# patient data 2 list
for(i in 1:length(patient_files)){
    DataList[[i]] = ROI2List(patient_files[i])
}
# health data 2 list
for (i in 1:length(health_files)){
    DataList[[length(patient_files)+i]] = ROI2List(health_files[i])
}

# Test, null hypothesis: the networks in the two groups are the same.
#load necessary libraries
library(RcppArmadillo)
library(Rcpp)
library(MTS)
library(Matrix)
library(MASS)
library(matrixcalc)
library(statmod)
sourceCpp('C:/Users/Clancy/Desktop/longitudinalFC-master/LongitInclude.cpp')
NPerm=5000
N1=length(patient_files)
N2=length(health_files)
CSFit=FCanalysis(datalist=DataList, N1=N1, N2=N2, Nperms=NPerm)
CSFit$T_global
permp((1-ecdf(CSFit$T_dist_global)(CSFit$T_global))*NPerm,NPerm,N1,N2)
Q = choose(length(interested_roi),2)
# multiple correlation
PLocal=rep(NA,Q)
for(q in 1:Q){
    PLocal[q]=permp((1-ecdf(CSFit$T_dist_local[,q])(CSFit$T_local[q]))*NPerm,NPerm,N1,N2)
}
Padj=round(p.adjust(PLocal,method="BH"),4)
Padj
