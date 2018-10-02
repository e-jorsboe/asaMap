args<-commandArgs(trailingOnly = T)

if(length(args)==0){
  
  print("Arguments have to be supplied: ")
  print("1. output file from asaMap")
  q()
}

fileName<-args[1]

fileName1<-tail(unlist(strsplit(fileName,"/")),1)

df<-read.table(fileName,as.is=TRUE,header=TRUE)

getP<-function(col1,col2,df){
  return(1-pchisq(-2*(col1-col2),df=df))
}

## add
if(ncol(df)==16){
  
  Ps<-cbind(df$Chromo,df$Position,df$nInd,df$f1,df$f2,getP(df$llh.M0.,df$llh.M1.,1),getP(df$llh.M1.,df$llh.M5.,2),
            getP(df$llh.M1.,df$llh.M2.,1),getP(df$llh.M1.,df$llh.M3.,1),getP(df$llh.M1.,df$llh.M4.,1),getP(df$llh.M2.,df$llh.M5.,1),
            getP(df$llh.M3.,df$llh.M5.,1),getP(df$llh.M4.,df$llh.M5.,1))
  
  colnames(Ps)<-c("Chromo","Position","nInd","f1","f2","M0vM1","M1vM5","M1vM2","M1vM3","M1vM4","M2vM5","M3vM5","M4vM5")
  
} else if(ncol(df)==23){
  
  Ps<-cbind(df$Chromo,df$Position,df$nInd,df$f1,df$f2,getP(df$llh.R0.,df$llh.R1.,1),getP(df$llh.R1.,df$llh.R7.,2),
            getP(df$llh.R1.,df$llh.R4.,1),getP(df$llh.R1.,df$llh.R5.,1),getP(df$llh.R1.,df$llh.R6.,1),getP(df$llh.R2.,df$llh.R6.,1),
            getP(df$llh.R3.,df$llh.R6.,1),getP(df$llh.R4.,df$llh.R7.,1),getP(df$llh.R5.,df$llh.R7.,1),getP(df$llh.R6.,df$llh.R7.,1))
  
  
  colnames(Ps)<-c("Chromo","Position","nInd","f1","f2","R0vR1","R1vR7","R1vR4","R1vR5","R1vR6","R2vR6","R3vR6","R4vR7","R5vR7","R6vR7")
  
} else{
  stop("Ouput file does not have expected dimensions!")
}
  
write.table(Ps,paste0(fileName1,".Pvalues"),col=T,row=F,qu=F)