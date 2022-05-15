# ==== Calculate BIC statistic for asReml-R objects  ==== #
#' @param x an asReml object from the call to asreml() 
#' @return integer of length equal to 1
#' \describe{
#' \item{x}{A numeric BIC value for the model fit in asReml}
#' }
library(asreml)

BIC_asREML<-function(x)  {log(x$nedf)*(length(x$vparameters)) - (2*x$loglik)
}