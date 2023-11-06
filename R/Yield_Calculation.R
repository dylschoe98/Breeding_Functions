#The purpose of the current script is to calculate maize grain yield from plot weight based on assuming:
#..the calculation assumes that the plot is 2 row with 20ft long rows and 30" row spacing but user-provided plot dimensions are allowed
#..assuming there are 56lbs of corn in a standard bushel 
#..for 2 row plot assume: 2 rows * 20ft long * 2.5ft wide = 100ft^2
#..43,560ft per acre

# ==== Calculate Neighbor Plant Height in a Field Trial ==== #
#' @param dt a data.table with a column for plot weight and grain moisture
#' @param Weight a vector of length 1 giving the name of the column that contains the plot weight information
#' @param Moisture a vector of length 1 giving the name of the column that contains the grain moisture information
#' @param Dimensions a vector of length 1 giving the plot area. Default is 100/43560 for assuming a 2 row plot of 100ft^2 
#' @return Data.table containing:
#' \describe{
#' \item{clean_dt}{ a data.table with all columns of the provided data data plus a column called 'GY' that calculates grain yield
#' }
#' 
library(data.table)
library(janitor)

#Take a date table with a column for  moisture and weight then return grain yield. Returns all intermediate steps in the corresponding calculation 
GrainYield<-function(clean_dt,Weight,Moisture,Dimensions=(100/43560)){
  plotDim<-Dimensions
  weight<-Weight
  moisture<-Moisture
    setnames(clean_dt,c(weight,moisture),c('weight','moisture'))
#1) Dry matter
  clean_dt$moist_account<-clean_dt$weight * (1 - (clean_dt$moisture * 10^-2))

#2) Dry matter to 15% moisture 
  clean_dt$dry_matter<- clean_dt$moist_account / 0.85 
  
#3) Bushels of corn harvested
  clean_dt$bushels <-  clean_dt$dry_matter / 56

#4) Bushels per acre ###  
  clean_dt$GY<-(clean_dt$bushels / plotDim)
return(clean_dt[,-c('moist_account','dry_matter','bushels')])}