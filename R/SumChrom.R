library(data.table)

SumChrom<-function(x){
  chr<-seq(1,10,1)    
  chr_new<-vector("list",length = length(chr))
  axisdf<-data.table(Chr=chr)
  mid.chr<-max.chr<-min.chr<-vector(length = length(chr))
  
  for (i in 1:length(chr)){
    chr_new[[i]]<-x[Chr==chr[i]]
    if (chr[i]==1){
      chr_new[[i]]$tot_pos<-x[Chr==chr[i]]$Pos
    }else{
      cum_chrs<-rbindlist(chr_new,fill = TRUE)
      cum_chrs[Chr==chr[i]]$tot_pos<-cum_chrs[Chr==chr[i]]$Pos + max(cum_chrs[Chr==chr[i-1]]$tot_pos)
      chr_new[[i]]<-cum_chrs[Chr==chr[i]]
    }
#Obtain the midpoint of the chromosome
    mid.chr[i]<-mean(chr_new[[i]][Chr==chr[i]]$tot_pos)
    
#Calculate min and max chromosome length
    chr_new[[i]][Chr==chr[i],ChromMin:=min(chr_new[[i]][Chr==chr[i]]$tot_pos)]
    chr_new[[i]][Chr==chr[i],ChromMax:=max(chr_new[[i]][Chr==chr[i]]$tot_pos)]
  }
  axisdf$mid.chr<-mid.chr
  return(list(rbindlist(chr_new),axisdf))
}