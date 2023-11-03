# ==== Calculate minor allele frequency (MAF) from a 0 & 1 marker matrix ==== #
#' @param x a genotype marker matrix of 0 & 1s where the first column corresponds to the name of the genotypes
#' \describe{
#' \item{maf_dt}{A data.table with two columns 'MAF', corresponding to the minor allele frequency of each marker and 'Marker' for the name of the genetic maker}
#' }

library(data.table)

MAF<-function(x){
#..convert to matrix
    colnames(x)[1]<-'Inbred'
    MarkerMatrix<-as.matrix(x[,c(-1)])
    rownames(MarkerMatrix)<-x$Inbred #giving the names of the individuals

#calculate MAF
maf_manual<-colSums(MarkerMatrix,na.rm = T) / nrow(MarkerMatrix) #going to ignore missing values
  maf_dt<-data.table(Marker=c(names(maf_manual)),MAF=c(maf_manual))
  maf_dt[MAF >0.5,MAF:=1-maf_dt[MAF >0.5]$MAF]
return(maf_dt)}

