#The purpose of the current script is to simualte DH gametes
# ==== Calculate Neighbor Plant Height in a Field Trial ==== #
#' @param dt a data.table with five columns: Genotype, Environment, Range, Row, and PH (values for plant height)... order and name of the columns does not matter. 
#' @return Data.table containing:
#' \describe{
#' \item{dtPH}{ a data.table with the first column corresponding to 'Genotype' for names of the genotypes provided,Row and Range from the original data, PHavg, and PHsum giving the average difference and sum of the cumulative difference between the two adjacent neighbors }
#' }
library(data.table)

NeighborPH<-function(x=dt){
  COL<-colnames(x)
    setnames(x,c(COL),c('Genotype','Environment','Range','Row','PH'))
  
    env<-unique(x$Environment)
#Loop over the ranges
  ListPH<-list()  
for (i in 1:length(env)){
    RA<-as.numeric(unique(x[Environment==env[i]]$Range))
  for (j in 1:length(RA)){
      RO<-as.numeric(unique(x[Environment==env[i] & Range==RA[j]]$Row))
  for (k in 1:length(RO)){
      temp<-x[Environment==env[i] & Range==RA[j]]
        temp$Row<-as.numeric(temp$Row)
        
#Subtract off neighboring distance
#..For the first plot in row
  if(RO[k]==min(RO)){
    subplot<-RO[k]+1
      if(!subplot %in% RO){DIF<-NA}else{
        DIF<-temp[Row==RO[k]]$PH - temp[Row==RO[k]+1]$PH
  }
}       
#..For the last plot in a row
  if(RO[k]==max(RO)){
    subplot<-RO[k]-1
      if(!subplot %in% RO){DIF<-NA}else{
        DIF<-temp[Row==RO[k]]$PH - temp[Row==RO[k]-1]$PH
  }
}       
#..For plots not the first or last        
if(RO[k] !=min(RO) | RO[k]!=max(RO)){
  subplot1<-RO[k]+1
    subplot2<-RO[k]-1
      DIF<-c(temp[Row==RO[k]]$PH - temp[Row==subplot1]$PH,temp[Row==RO[k]]$PH - temp[Row==subplot2]$PH)
}       
        
#List out Values
  tempdt<-cbind(temp[Row==RO[k]][,-c('PH')],data.table(PHavg=c(mean(DIF,na.rm = T)),PHsum=c(sum(DIF,na.rm=T)),PHmax=c(max(DIF,na.rm=T))))
ListPH[[paste0(i,'_',j,'_',k)]]<-tempdt
    }   
  }
}
return(dtPH=rbindlist(ListPH))}
