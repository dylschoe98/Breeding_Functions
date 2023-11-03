#The purpose of the current script is to calculate rep-to-rep correlations within a plant breeding field trial.

#' @param x a data.table with columns at minimum for genotype names, replictation, and the phenotypic values
#' @param pedigree a vector of length one corresponding to the name of the column for genotypes
#' @param Rep a vector of length one corresponding to the name of the column for replicates
#' @param method a vector of length one corresponding to either perason, kendall and spearman to provide the method of correlation
#' @param use a vector of length one corresponding to eithereverything, all.obs, complete.obs, na.or.complete, pairwise.complete.obs for the function cor. Default is pairwise.complete.obs

#' @return a data.table
#' \describe{
#' \item{corDT}{A data.table with a column Desicription and Cor where Cor provides the correlation between the replicates}
#' }

library(data.table)

ReptoRep<-function(x,pedigree, Rep,method=c("pearson", "kendall", "spearman"),use='pairwise.complete.obs'){
  PED<-pedigree
  REP<-Rep
  METH<-method
  USE<-use
  
#Ensure proper specifications     
if(!METH %in% c("pearson", "kendall", "spearman")){
    stop('method must be either pearson, kendall, or spearman')
}    
  
if(!USE %in% c("everything", "all.obs", "complete.obs", "na.or.complete", "pairwise.complete.obs")){
    stop('use must be either everything, all.obs, complete.obs, na.or.complete, pairwise.complete.obs')
}    
  
#Dcast and correlate    
yDT<-dcast.data.table(x,as.formula(paste0(PED,'~',REP)),value.var = 'value')
  colnames(yDT)<-c('pedigree','Rep1','Rep2')
    corDT<-data.table(Description=c('Correlation'),Cor=cor(x = yDT$Rep1,yDT$Rep2,method = METH,use=USE))
return(corDT)}