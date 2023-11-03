# ==== Perform K-fold Cross Validation (CV) for GP/GS Models using ridge regression  ==== #
#' @param dt data.table with one response for each genotype. The column names should be 'genotype' & 'value'
#' @param G_List a list of relationship matricies of NxN dimensions where N is the number of genotypes and M is the number of markers. Rownames need to be provided to the matrix with the genotype 
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
CV<-function(dt,G_List,k=5,Repeat=20){
  ObjectName <- function(x) deparse(substitute(x)) #this is a function used to obtain the name of the object that is the matrix such that it can be save in a list
  
ResultsList<-list()
REPS<-Repeat
k

for (r in 1:REPS){ #for each rep, use the k-1 groups as a trainings set to predict the kth group
  
#Randomly shuffle the data
data_shuffle<-dt[sample(1:nrow(dt)),]
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
    for (g in 1:length(G_List)){
          G<-G_List[[g]] #Subset out the first matrix  
          
#Sample 100% -k% for training      
    TrainSet<-data_shuffle[Group !=GROUPS[i]]
      TestSet<-data_shuffle[!genotype %in% TrainSet$genotype]
    CompleteSet<-rbind(TrainSet,TestSet) 

#subset out the phenotypes 
  y<-CompleteSet$value
  names(y)<-CompleteSet$genotype
    y[names(y) %in% TestSet$genotype]<-NA
#...add NAs to the y values     
      y<-y[order(factor(names(y), levels = c(rownames(G))))]
          
# Run rrBLUP
  ans_real_BV <- mixed.solve(y,K=G) 
    print(paste('Rep:',r,'Group',GROUPS[i],'Matrix Number',g))
          
# Summarize the values 
bv_real<-ans_real_BV$u + as.numeric(ans_real_BV$beta)
          
#..Obtain missing values
y_gebv<-bv_real[names(bv_real) %in% TestSet$genotype]
  y_gebv<-y_gebv[sort(names(y_gebv))]
#...blues
     y_blues<-TestSet$value
      names(y_blues)<-TestSet$genotype
        y_blues<-y_blues[sort(names(y_blues))]
ResultsList[[paste0(g,'_',GROUPS[i],'_',r)]]<-data.table(Group=c(GROUPS[i]),Cor=c(cor(y_blues,y_gebv)),Matrix=c(paste0(ObjectName(G),'_',g)))
    }
  }
}
Pred_DT<-rbindlist(ResultsList)[,lapply(.SD,mean),by=.(Matrix,Group),.SDcols=c(2)]
return(Pred_DT)}

