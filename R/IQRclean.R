# ==== Calculate Upper and Lower Limits based on an IQR Test  ==== #
#' @param x a list of vector of values to identify the upper and lower limits 
#' @param threshold the value to multiple the IQR by to establish upper and lower limits for removal
#' @return a list that each contains a vector of length 2
#' \describe{
#' \item{quantiles_removed}{a vector of length 2 with upper and lower limits}
#' }
IQRclean<-function(x,threshold=2.5){
  Thresh<-threshold
  quantiles_removed<-list()
for (i in 1:length(x)){
    quantiles_removed[[i]]<-c(quantile(x[[i]],0.75,na.rm=T)+Thresh*IQR(x[[i]],na.rm=T),
                           quantile(x[[i]],0.25,na.rm=T)-Thresh*IQR(x[[i]],na.rm=T))
}  
return(quantiles_removed)}