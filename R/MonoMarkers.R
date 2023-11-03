# ==== Identify genetic markers that are monomorphic from a 0 & 1 marker matrix ==== #
#' @param x genotype matrix of 0 and 1s with the first column being the name of the genotypes
#' \describe{
#' \item{MonoMarkers}{a data.table with one column name, 'Markers,' corresponding to the name of the markers that are monomorphic}
#' }

library(data.table)
MonoMarkers<-function(x){
  Markers<-colnames(x[,-c(1)])
#..identify the SNPs that are Monomorphic in the parents
  mean_cols<-colMeans(x[, ..Markers],na.rm = T)
    monoFound<-mean_cols[mean_cols ==1]
      MonoMarkersDT<-data.table(Markers=c(names(monoFound))) 
return(MonoMarkersDT)}
