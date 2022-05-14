# ==== Perform K-fold Cross Validation (CV) for GP/GS Models using ridge regression  ==== #
#' @param dt data.table with one response for each genotype. The column names should be 'genotype' & 'value'
#' @param G a marker matrix of NxM dimensions where N is the number of genotypes and M is the number of markers. Rownames need to be provided to the matrix and NO missing values can be present in the marker matrix or rrBLUP will throw and error 
#' @param k a numeric value corresponding to the folds of CV to be complete. Default is 5 fold CV.
#' @param Repeat a numeric value corresponding to the number of times the CV is repeated per group. Default is 20 repetitions per fold.
#' @return Data.table containing
#' \describe{
#' \item{Pred_DT}{A data.table with a correlation for each 'k' group. Correlations are between the estimated true value from the Model and the BLUEs/BLUPs}
#' }
library(data.table)
library(rrBLUP)
library(dplyr)

#1) Perform k fold CV 
CV<-function(dt,G,k=5,Repeat=20){

ResultsList<-GroupList<-list()
REPS<-Repeat
k

#Randomly shuffle the data
data_shuffle<-pheno_sub[sample(1:nrow(dt)),]
  n<-nrow(data_shuffle)   
    group_assign<-rep(letters[1:k],each=n/k)
      n_diff<- n-length(group_assign)
      
#Ensure equal number of letters
#..more observations then assignments
  if(n_diff >0){group_assign<-c(group_assign,(letters[1:abs(n_diff)]))}
#..more assignments than observations
  if(n_diff <0){group_assign<-group_assign[1:(length(group_assign)-n_diff)]}
      
data_shuffle$Group<-group_assign
GROUPS<-unique(data_shuffle$Group)

for (i in 1:length(GROUPS)){  
  for (r in 1:REPS){
    
#Sample 100% -k% for training      
    TrainSet<-sample_n(data_shuffle,size = nrow(data_shuffle[Group !=GROUPS[i]]))
      TestSet<-data_shuffle[!genotype %in% TrainSet$genotype]
    CompleteSet<-rbind(TrainSet,TestSet) 

#subset out the phenotypes 
  y<-CompleteSet$value
  names(y)<-CompleteSet$genotype
    y[names(y) %in% TestSet$genotype]<-NA
#...add NAs to the y values     
      y<-y[order(factor(names(y), levels = c(rownames(G))))]
          
# Run rrBLUP
  ans_real_BV <- mixed.solve(y,K=A.mat(G)) 
    print(paste('Rep:',r,'Group',GROUPS[i]))
          
# Summarize the values 
bv_real<-ans_real_BV$u + as.numeric(ans_real_BV$beta)
          
#..Obtain missing values
y_gebv<-bv_real[names(bv_real) %in% TestSet$genotype]
  y_gebv<-y_gebv[sort(names(y_gebv))]
#...blues
     y_blues<-TestSet$value
      names(y_blues)<-TestSet$genotype
        y_blues<-y_blues[sort(names(y_blues))]
ResultsList[[paste0(GROUPS[i],'_',r)]]<-data.table(Group=c(GROUPS[i]),Cor=c(cor(y_blues,y_gebv)))      
  }
GroupList[[paste0(GROUPS[i],'_',r)]]<-rbindlist(ResultsList)[,lapply(.SD,mean),by=.(Group),.SDcols=c(3)]
}
return(Pred_DT=rbindlist(GroupList))}

