# ==== Calculate BLUEs or BLUPs, Variance Components and Heritability  ==== #
#' @param Mod an asreml or mmer object from the Sommer R package 
#' @param Mod.Name a vector for for the name desired to be called to the model 
#' @return A data.table containing
#' \describe{
#' \item{data}{A data.table with the variance components for the random effects of the linear model. The column 'Proportion' refers to the prroportion of phenotypic variance explained by each random component}

library(data.table)

VC.datatable<-function(Mod,Mod.Name){
  Mod.Name
  
  data<-as.data.table(summary(Mod)$varcomp)
  Mod.Check<-grep("component",names(data)) #to distinguish if asREML or if mmmer was used
  
  if(length(Mod.Check) != 0){
    data2<-summary(Mod)$varcomp
    data$Components<-rownames(data2)
    colnames(data)[1]<-"VarComp"
    data<-data[,-(2:5)]
    data[Components=="units!R",Components:="Residual"]
  }else{
    ran.term<-unlist(Mod$terms$random)
    res.term<-unlist(Mod$terms$rcov)
    mod.terms<-append(ran.term,res.term)
    data$Components<-mod.terms
    data<-data[,-(2:4)]
    data[Components=="units",Components:="Residual"]
  }  
  data$Mod.Name<-Mod.Name
  data$Proportion<-data$VarComp / sum(data$VarComp)
return(data)  
}
