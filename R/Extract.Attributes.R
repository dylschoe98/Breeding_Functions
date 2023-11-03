#The purpose of the current script is to attract attributes from B73v4 and B73v5 gene attotation file
#' @param gtf a data.table corresponding to the attribute file
#' @param attribute a vector of length one corresponding to the name of the column for subsetting out
#' @param col.name a vector of length one corresponding to the name for the new column 
#' @param Version a vector of length one either v4 or v5 depending on the version of the attribute file
#' @return a list
#' \describe{
#' \item{formated.genes}{A data.table with with a new column that is provided by col.name and corresponds to the attributes of interest}
#' }

Extract.Attribute<-function(gtf, attribute,col.name,Version){ #Version can take v4, v5, or description
  library(stringr)
  attribute
  col.name
  gtf<-as.data.table(gtf)
  sub.gtf<-gtf[attributes %like% attribute]
  
if(Version=="v4"){sub.gene<-as.data.table(str_match(sub.gtf$attributes, paste(attribute,"=(.*?);",sep = "")))}
  if(Version=="v5"){sub.gene<-as.data.table(str_match(sub.gtf$attributes, paste(attribute," (.*?);",sep = "")))}  
    if(Version=="description"){sub.gene<-as.data.table(str_match(sub.gtf$attributes, paste(attribute,"=(.*?);",sep = "")))}
#..add back in the gene id as a column
gene.id.Ver<-cbind(sub.gtf,sub.gene[,2])
  colnames(gene.id.Ver)[10]<-"gene"
    ck<-as.data.frame(sapply(gene.id.Ver, function(x) gsub("\"", "", gene.id.Ver$gene)))
formated.genes<-cbind(gene.id.Ver[,1:9],ck[,10])
  colnames(formated.genes)[10]<-col.name
return(formated.genes)}   