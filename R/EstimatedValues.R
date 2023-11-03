# ==== Calculate BLUEs or BLUPs, Variance Components and Heritability  ==== #
#' @param dt data.table with a column corresponding to at least each term in the desired mixed model. Each term in the model must properly be specified as a factor. Also a column for the phenotype identifier must be provided as 'variable'
#' @param Pheno a vector of phenotype type names for the mixed model
#' @param Geno the name of the column that corresponds to the genotype term in the linear model
#' @param Random a vector for the desired random effects of the model
#' @param Fixed a vector for the desired fixed effects of the model
#' @param Residual an optional vector for the residual structure. If not vector is provided the default is independent and identically distributed resiudal structure: '~id(units)' 
#' @param Weight the name of the column if the residuals should be weighted. The default is FALSE
#' @param ReturnWeights a TRUE or FALSE statement telling whether the weight of each data point calculated from 1/sqrt(SE) should be returned as described by Smith et al., 2001
#' @return List containing
#' \describe{
#' \item{ExpVal_dt}{A data.table with a BLUE or BLUP corresponding to each genotype per pheno}
#' \item{VC_dt}{A data.table with the variance components of the random terms in the model per pheno}
#' \item{H2_dt}{A data.table with the heritability calculated using the Cullis Method (Cullis et al., 2006) per pheno}
#' }
library(data.table)
library(asreml)

EstimatedValues<-function(dt,Pheno,Geno, Random,Fixed,Residual='~id(units)',Weight=FALSE,ReturnWeights=FALSE){

source("~/Breeding_Functions/R/VC.data.table.R") #must call the internal function first
  
Weight<-Weight
S2.Weight<-Weight
Geno<-Geno
Ran<-Random
Fix<-Fixed
Resid<-Residual
phenos<-Pheno
ExpVal_List<-VC_List<-H2_List<-list()

for (i in 1:length(phenos)){
  pheno_sub<-dt[variable==phenos[i]]
  pheno_copy<-copy(pheno_sub)    
  
#Call to Asreml for Non-Weighted Analysis or Weighted Analysis
if(Weight==FALSE){
  Mod<-asreml(fixed = as.formula(paste("value~",Fix)), 
                  random=as.formula(paste0("~",Ran)),residual =  as.formula(paste("~",Resid)),na.action = na.method(x=c("include")),data = pheno_sub)
  if(!Geno %like% Ran){
    FixSub<-gsub(paste0(Geno),'1',Fix)
#Make a genotype a random term for heritability
  ModRand<-asreml(fixed = as.formula(paste("value~",FixSub)), 
              random=as.formula(paste0("~",Geno,'+',Ran)),residual = as.formula(paste("~",Resid)),na.action = na.method(x=c("include")),data = pheno_copy)
  
#Random Predicted and variance components
PredRand<-predict.asreml(ModRand,classify = Geno, present = Geno)
  vc<-VC.datatable(ModRand,Mod.Name = phenos[i])
H2_List[[i]]<-data.table(H2_Method=c('Cullis'),Pheno=c(phenos[i]),H2=c((1 - ((PredRand$avsed**2) / (2 * vc[Components==Geno]$VarComp)))* 100))
  } #end the second if statement
}else{
  Mod<-asreml(fixed = as.formula(paste("value~",Fixed)), 
              random=as.formula(paste0("~",Random)),residual = as.formula(paste("~",Resid)),
                family=asr_gaussian(dispersion=1),weights = as.formula(paste0(S2.Weight)),
                    na.action = na.method(x=c("include")),data = pheno_sub)
  if(!Geno %like% Ran){
      Fix<-gsub(paste0(Geno),'1',Fix)
#Make genotype a random term for heritability
  ModRand<-asreml(fixed = as.formula(paste("value~",Fix)), 
                random=as.formula(paste0("~",Geno,'+',Ran)),residual = as.formula(paste("~",Resid)),
                    family=asr_gaussian(dispersion=1),weights = as.formula(paste0(S2.Weight)),
                       na.action = na.method(x=c("include")),data = pheno_copy)
#Random Predicted and variance components
PredRand<-predict.asreml(ModRand,classify = Geno, present = Geno)
  vc<-VC.datatable(ModRand,Mod.Name = phenos[i])
H2_List[[i]]<-data.table(H2_Method=c('Cullis'),Pheno=c(phenos[i]),H2=c((1 - ((PredRand$avsed**2) / (2 * vc[Components==Geno]$VarComp)))* 100))
  }
} #end the else statement

#Return the BLUEs or BLUPs
  Pred<-predict.asreml(Mod,classify = Geno, present = Geno)
  ExpVal<-as.data.table(Pred$pvals)
  ExpVal$Pheno<-phenos[i]
  
#Return the BLUEs or BLUPs
  if(Geno %in% Ran){
#..Heritability 
    H2_List[[i]]<-data.table(H2_Method=c('Cullis'),Pheno=c(phenos[i]),H2=c(1 - ((Pred$avsed**2) / (2 * vc[Components==Geno]$VarComp))* 100))
  }
  
#..For weights
  if(ReturnWeights==TRUE){
      ExpVal$S2.weight<-1 / (sqrt(ExpVal$std.error))}
ExpVal_List[[i]]<-ExpVal

#Variance Components
VC_List[[i]]<-VC.datatable(Mod,Mod.Name = phenos[i])
}
return(list(ExpVal_dt=rbindlist(ExpVal_List),VC_dt=rbindlist(VC_List),H2_dt=rbindlist(H2_List)))}

