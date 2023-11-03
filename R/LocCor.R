#The purpose of the current script is to calculate correlations between environments within a plant breeding MET.
#' @param x a data.table with columns at minimum for genotype names, replictation, and the phenotypic values. Phenotypic values must be specified with the column name 'value'.
#' @param pedigree a vector of length one corresponding to the name of the column for genotypes
#' @param Environment a vector of length one corresponding to the name of the column for environments
#' @param Phenotypes a vector of length one corresponding to the name of the column for phenotypes
#' @param method a vector of length one corresponding to either perason, kendall and spearman to provide the method of correlation
#' @param use a vector of length one corresponding to eithereverything, all.obs, complete.obs, na.or.complete, pairwise.complete.obs for the function cor. Default is pairwise.complete.obs

#' @return a list
#' \describe{
#' \item{corDT}{A data.table with all pair-wise correlations between environments per phenotype}
#' \item{corList}{A list of correlation correlation matricies}
#' \item{pList}{A list of p-value matricies}
#' }

library(data.table)
library(ggcorrplot)

LocCor<-function(x,pedigree, Environment,Phenotypes,method=c("pearson", "kendall", "spearman"),use='pairwise.complete.obs'){
  PED<-pedigree
  ENV<-Environment
  PHE<-Phenotypes
  METH<-method
  USE<-use
  
#Ensure proper specifications     
if(!METH %in% c("pearson", "kendall", "spearman")){
    stop('method must be either pearson, kendall, or spearman')
}    
  
if(!USE %in% c("everything", "all.obs", "complete.obs", "na.or.complete", "pairwise.complete.obs")){
    stop('use must be either everything, all.obs, complete.obs, na.or.complete, pairwise.complete.obs')
}    
 
COLS<-colnames(x)[colnames(x)=='value']
if(length(COLS) !=1){
    stop('Phenotyic values must be specified by the column name value.')
  }    
  
#Dcast and correlate    
yDT<-dcast.data.table(x,as.formula(paste0(PED,'+',PHE,'~',ENV)),value.var = 'value')
  setnames(yDT,c(PED,PHE),c('pedigree','phenotype'))

#Correlations
phenos<-unique(as.character(yDT$phenotype))  
CorList<-P<-CorPairwise<-list()
for (i in 1:length(phenos)){
  
P[[i]]<-cor_pmat(yDT[phenotype==phenos[i]][,-c('pedigree','phenotype')])
  CorList[[i]]<-cor(yDT[phenotype==phenos[i]][,-c('pedigree','phenotype')],use = use,method = METH)
#..pairwise combinations between environments
  corMelt<-as.data.table(reshape2::melt(cor(yDT[phenotype==phenos[i]][,-c('pedigree','phenotype')],use = use,method = METH)))
    corMelt$EnvPairs<-paste0(corMelt$Var1,'_',corMelt$Var2)
      setnames(corMelt,c('value'),c('r'))
#..add in the pvalues
    pvalmelt<-as.data.table(reshape2::melt(cor_pmat(yDT[,-c('pedigree','phenotype')])))
      setnames(pvalmelt,c('value'),c('Pvalue'))
      corMelt<-cbind(corMelt,pvalmelt[,c('Pvalue')])
#..do not include duplicate    
    corMelt<-corMelt[,.SD[which.max(r)],by=.(r)]
      corMelt$Phenotype<-phenos[i]
        CorPairwise[[i]]<-corMelt
  }    
  corDT<-rbindlist(CorPairwise)
    names(CorList)<-phenos
return(list(CorList,P,corDT))}
