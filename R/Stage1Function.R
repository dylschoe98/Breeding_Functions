# ==== Stage 1 Analysis ==== #
#' @param dt.tidy dataframe with a environment, variable (phenotype names), rep,value (phenotype response), pedigree, tester, range, and pass 
#' @param dt.covar dataframe with columns for environment, rep, and all phenotype names as column names that will be used as a covariate along with a column for range and pass.
#' @param covars vector of model covariates 
#' @param no_covar vector of phenotypes that will not be modeled via a covariate 
#' @return List containing
#' \describe{
#' \item{dt}{A seven column dataframe with  'pedigree', 'predicted.value', 'std.error', 'status', 'pheno','S2.weight','environment'}
#' \item{h2}{A dataframe with a count of the number of times and inbred made the top of the list}
#' }

library(data.table)
library(asreml)
library(dplyr)

# asreml.options(ai.sing=FALSE)

#Function to remove NAs from the Variance-Covariance matrix and obtain diagonal weights for the Stage 2 analysis  
S2.Weights_asreml<-function(mat){
  NA_list<-apply(mat,MARGIN = 2,function(x){
    col.x<-x[!(is.na(x))]
    ck.dt<-as.data.table(col.x)
    return(ck.dt)})
  
  No.NA_list<-list()
  for(i in 1:length(NA_list)){
    if(nrow(NA_list[[i]]) >0)
      No.NA_list[[i]]<-NA_list[i]
  }
  ck2<-bind_cols(No.NA_list)
  mat_new<-as.matrix(ck2) 
#..get the weights 
  square_noNA<-solve(mat_new)  
  S2.Weights<-diag(square_noNA)
return(S2.Weights)}


#Practice
# dt.tidy<-copy(cleanDT[environment=='INH1_2021' &tester=='PHK76' & variable=='RL'])
# DF<-30
# SpatCor<-0.60
# no_covar<-NoCovar
# dt.covar<-copy(covar.Clean[environment=='INH1_2021' &tester=='PHK76'])
# TraceASREML=FALSE
# covars<-colnames(covar.Clean)[8:14]

Stage1Analysis<-function(dt.tidy,dt.covar,covars=NULL,no_covar=NULL,SpatCor=0.60,DF=30,TraceASREML=FALSE){
  asreml.options(gammaPar=TRUE)
  
if(!is.numeric(DF)){
  stop('The degrees of freedom must be numeric.')
}
    
asreml.options(trace=TraceASREML)
  BIC<-function(x)  log(x$nedf)*(length(x$vparameters)) - (2*x$loglik) #from Jim Holland as published in Rogers et al. (2021) for calculating Heritability
  
locs<-unique(as.character(dt.tidy$environment))
  VarCompList<-WaldTestList<-CovList<-BLUEsList<-H2ResidualList<-H2CovariateList<-VcovList <-skip.list<-ModelBIC<-list()
for (j in 1:length(locs)){
    
loc_sub<-dt.tidy[environment==locs[j]]
testers<-unique(as.character(loc_sub$tester))
  for (t in 1:length(testers)){
        test.sub<-loc_sub[tester==testers[t]]
    
#phenotypes of interest          
pheno<-as.character(unique(factor(test.sub$variable))) #pheno type list
          
#Qualtiy control check ensuring at least a set number of DF and variane 
  DF.check<-test.sub[, .N, by=.(variable)]

#Format
test.sub$variable<-factor(test.sub$variable)
test.sub$pedigree<-factor(test.sub$pedigree)
    for (i in 1:length(pheno)){
      
#...subset out the phenotypes
    sub.pheno<-test.sub[variable==pheno[i]]
      sub.pheno$rep<-as.integer(sub.pheno$rep)
      
#QC on the data 
      
#...if total observations < DF, skip
  if(nrow(sub.pheno[variable==pheno[i] & !is.na(value)]) < DF) {skip.list[[paste0(locs[j],'_', testers[t],'_', pheno[i])]] = data.table(pheno=c(pheno[i]), 
                                                                                                 tester=paste(testers[t]), environment=paste(locs[j]),Resason=c(paste0("<",DF,"observations")))
        print(paste0('Not enough DF for ',locs[j],'_', testers[t],'_', pheno[i], ' moving to the next phenotype for ',locs[j]))
      next} 
#...if values are all NA or have a variance of 0
  if(var(sub.pheno[variable==pheno[i]]$value,na.rm = T)==0 |is.na(var(sub.pheno[variable==pheno[i]]$value,na.rm = T))){ 
    skip.list[[paste0(locs[j],'_', testers[t],'_', pheno[i])]] = data.table(pheno=c(pheno[i]), tester=paste(testers[t]), environment=paste(locs[j]),Resason=c("Variance equal to 0"))
      print(paste0('0 Variance for ',locs[j],'_', testers[t],'_', pheno[i], ' moving to the next phenotype for ',locs[j]))
    next}    
      
#Format the data        
sub.pheno$rep<-as.integer(sub.pheno$rep)

#1) Covariate Tests ####         

  
#A) Test a Model without a covariate #### 
Mod.No_Covar<-tryCatch(asreml(fixed = value ~ pedigree, random= ~rep, data = sub.pheno[variable==pheno[i]]),
                       error=function(x) {print("Failed") 
                          return(NA)})

#Check if the model is NA, if NA restart the loop, use the length for this step 
  Model_Check<-copy(Mod.No_Covar)
if(is.na(Model_Check$loglik)){
  skip.list[[paste0(locs[j],'_', testers[t],'_', pheno[i])]] = data.table(pheno=c(pheno[i]), tester=paste(testers[t]),
                                                environment=paste(locs[j]),Resason=c("Model failed"))
    print("Model failed, restart loop")
      next
}

#Test each of the covariates if the phenotype was not in the no_covar vector        
  if(!(pheno[i]) %in% no_covar & !is.null(no_covar)){
        CovarSub<-dt.covar[environment==locs[j] & tester==testers[t]]
#..Add the covariates to the data 
    prep_w10004<-merge(sub.pheno,CovarSub, by.x = c("environment","pedigree","tester","rep",'range','pass'),by.y = c("environment","pedigree","tester","rep",'range','pass'))
      setnames(prep_w10004,c('value'),c('pheno_value')) #For the melt to work will need to fix the phenotype name
        w10004_sub_melt<-melt(prep_w10004,id.vars = c("environment","pedigree","rep",'range','pass','pheno_value'),measure.vars = c(colnames(prep_w10004[, ..covars])))
        
#Format the data         
  w10004_sub_melt$environment<-as.factor(w10004_sub_melt$environment)
    w10004_sub_melt$pedigree<-as.factor(w10004_sub_melt$pedigree)
#..Number of reps    
    REPS<-length(unique(w10004_sub_melt$rep))
    
#..B) Significance of covariates ####
lt_covar<-unique(as.character(w10004_sub_melt$variable))

Wald.Test.CovarianceList<-VarComp.Covariance_List<-vector("list", length = length(covars))
names(Wald.Test.CovarianceList)<-names(VarComp.Covariance_List)<-covars

#Test each of the covariates    
  for (k in 1:length(lt_covar)){
#..Model as fixed    
    Mod.Covar<-tryCatch(asreml(fixed = pheno_value~ pedigree+value, random= ~rep,na.action = na.method(x=c("include")),
                      data = w10004_sub_melt[variable==lt_covar[k]]),error=function(x)return(NA))
      if(is.na(Mod.Covar$loglik)){next} #restart the loop if the model fails
          Wald.Test.Covariance<-wald(Mod.Covar)
#...Save the Fixed effect Results      
      Wald.Test.Covariance_DT<-as.data.table(as.data.frame(Wald.Test.Covariance),keep.rownames = "id")
        Wald.Test.Covariance_DT[id=="value",id:=paste(lt_covar[k])]
          Wald.Test.Covariance_DT$Covariate<-lt_covar[k]
      Wald.Test.Covariance_DT$BIC<-BIC(Mod.Covar)
        Wald.Test.CovarianceList[[k]]<-Wald.Test.Covariance_DT
        
#..B) Heritability of model covariate ####
        
#..Model as random
    Mod.Covar_ran<-tryCatch(asreml(fixed = pheno_value~ value, random= ~pedigree+rep,na.action = na.method(x=c("include")), 
                    data = w10004_sub_melt[variable==lt_covar[k]]),error=function(x)return(NA))
        if(is.na(Mod.Covar_ran$loglik)){next} #restart the loop if the model fails s
      pred_ran_covar<-predict.asreml(Mod.Covar_ran,classify = "pedigree", present = "pedigree")
#...Save the variance components          
    VarComp.Covariance_DT<-as.data.table(summary(Mod.Covar_ran)$varcomp,keep.rownames="id")
      VarComp.Covariance_DT[id=="value",id:=paste(lt_covar[k])]
        VarComp.Covariance_DT$Covariate<-lt_covar[k]
    VarComp.Covariance_DT$H2_entry<-VarComp.Covariance_DT[id=="pedigree"]$component / (VarComp.Covariance_DT[id=="pedigree"]$component + (VarComp.Covariance_DT[id=="units!R"]$component/REPS)) *100
      VarComp.Covariance_DT$H2_cullis<-(1 - ((pred_ran_covar$avsed**2) / (2 * VarComp.Covariance_DT[id=="pedigree"]$component)))* 100
        VarComp.Covariance_DT$BIC<-BIC(Mod.Covar_ran)
          VarComp.Covariance_List[[k]]<-VarComp.Covariance_DT
rm(Mod.Covar_ran)        
  }
} #ending the covariate loop

#..C) Wald Test for the significance of Covariates ####    
if(!is.na(Model_Check$loglik)){  
  Wald.Test.No_Covariance<-wald(Mod.No_Covar)
    Wald.Test.No_Covariance_DT<-as.data.table(as.data.frame(Wald.Test.No_Covariance),keep.rownames = "id")
      Wald.Test.No_Covariance_DT$Covariate<-"None"
  Wald.Test.No_Covariance_DT$BIC<-BIC(Mod.No_Covar)
    rm(Mod.No_Covar)
#rbind with the covariates
    if(!(pheno[i]) %in% no_covar & !is.null(no_covar)){   
        Wald.Test_DT<-rbind(rbindlist(Wald.Test.CovarianceList),Wald.Test.No_Covariance_DT)
    }else{Wald.Test_DT<-Wald.Test.No_Covariance_DT}  
Wald.Test_DT[,`:=`(pheno=paste(pheno[i]), tester=paste(testers[t]), environment=min(locs[j]))]
      }
  
#..D) Test a Random model without a covariate ####
    Mod.NoCovar_ran<-asreml(fixed = value ~ 1, random= ~pedigree+rep,data = sub.pheno[variable==pheno[i]])
      pred_ran<-predict.asreml(Mod.NoCovar_ran,classify = "pedigree", present = "pedigree")
#..Number of reps
    REPS<-length(unique(sub.pheno$rep))
    
#Save the Results       
VarComp.No.Covariance_DT<-as.data.table(summary(Mod.NoCovar_ran)$varcomp,keep.rownames="id")
  VarComp.No.Covariance_DT$Covariate<-"None"
    VarComp.No.Covariance_DT$H2_entry<-(VarComp.No.Covariance_DT[id=="pedigree"]$component / (VarComp.No.Covariance_DT[id=="pedigree"]$component + (VarComp.No.Covariance_DT[id=="units!R"]$component)/REPS)) *100
VarComp.No.Covariance_DT$H2_cullis<-(1 - ((pred_ran$avsed**2) / (2 * VarComp.No.Covariance_DT[id=="pedigree"]$component)))* 100
  VarComp.No.Covariance_DT$BIC<-BIC(Mod.NoCovar_ran)
        
# Combine  the heritability by covariate 
#..Scenairos where the term is not in the no_covar object        
  if(!(pheno[i]) %in% no_covar & !is.null(no_covar)){ 
#...for scenaiors where no covariate is provided
      if(is.null(covars)==TRUE){Her.loc<-rbind(VarComp.No.Covariance_DT) 
        }else{ 
#...for scenaiors where covariates are provided
      Her.loc<-rbind(rbindlist(VarComp.Covariance_List),VarComp.No.Covariance_DT)}
#...for scenaiors where phenotype is in the no_covar object
  }else{Her.loc<-VarComp.No.Covariance_DT}  
  
#Identify the highest heritability and model usign BIC  
Her.loc$Best_BIC<-Her.loc[, .SD[which.min(BIC)]]$Covariate
  Her.loc$Best_H2<-Her.loc[, .SD[which.max(H2_entry)]]$Covariate
    Her.loc[, `:=`(environment=paste(locs[j]),tester=paste(testers[t]),pheno=paste(pheno[i]))]
    
#2) Spatial Analysis using range and pass information ####
          
#lists to store the results of the spatial analysis
  Spatial.results<-H2_tester<-VC.tester<-wald_stat.tester<-vector("list",length = length(testers)) #list to store the results
    names(VC.tester)<-names(H2_tester)<-names(wald_stat.tester)<-names(Spatial.results)<-testers
        
#..A) Significance of the Covariate for the spatial model ####
#guarantees the covariate is significant at an alpha value of 0.05 and if multiple significant covariate are detect the covariate with the lowest BIC is selected
  if(!(pheno[i]) %in% no_covar & !is.null(no_covar)){
      if(Wald.Test_DT[id %in% lt_covar,.SD[which.min(`Pr(Chisq)`)]]$`Pr(Chisq)` >= 0.05){sign.covar<-Wald.Test_DT[Covariate=="None"]
      }else{
        sign.covar<-Wald.Test_DT[`Pr(Chisq)`< 0.05 & id %in% covars]}
  }else{sign.covar<-Wald.Test_DT}
    sig.Her.loc<-Her.loc[Covariate %in% sign.covar$Covariate]
      Covariate<-sig.Her.loc[,.SD[which.max(H2_cullis)]]$Covariate  
    
#Spatial analysis for locations where the range and pass is provided   
  if(!(all(is.na(sub.pheno$range)) | all(is.na(sub.pheno$pass)))){
            
#....get all possible range and pass combinations
      R<-seq(min(sub.pheno$range), max(sub.pheno$range))
        P<-seq(min(sub.pheno$pass), max(sub.pheno$pass))
          r_p<-as.data.table(expand.grid(R,P))
      colnames(r_p)<-c("range","pass")
        r_p$environment<-unique(sub.pheno$environment)
          r_p$RP<-paste(r_p$range,"_",r_p$pass, sep = "")
          
#subset out the necessary column names
  if(!is.null(dt.covar)){
  CovarSub_Spat<-dt.covar[environment==locs[j] & tester==testers[t]]}
    if(!(pheno[i])%in% no_covar & Covariate !='None'){ #Such that if a covariate is not needed the loop moves to the else statement 
      
#..Add the covariates to the data 
  sub_clean<-merge(sub.pheno,r_p,by.x=c("environment","range","pass"), by.y = c("environment","range","pass"),all = T)
  sub_clean$pass<-as.integer(sub_clean$pass) 
    setnames(sub_clean,c('value'),c('pheno_value')) #Change the covariate name to value 
#...merge the complete set of range and columns with the covariates   
  COLS<-c("environment","pedigree","tester","rep",'range','pass',Covariate)
      merged.dt<-merge(sub_clean,CovarSub[, ..COLS], by.x = c("environment","pedigree","tester","rep",'range','pass'),
                         by.y = c("environment","pedigree","tester","rep",'range','pass'),all = T)
      setnames(merged.dt,c(Covariate),c('value')) #Change the covariate name to value 
    }else{
        merged.dt<-merge(sub.pheno,r_p,by.x=c("environment","range","pass"), by.y = c("environment","range","pass"),all = T)
      }      
          
#Order the data
  merged.dt<-merged.dt[order(merged.dt$range,merged.dt$pass,decreasing = F)]
  
#add the correct cases for asREML
merged.dt$environment<-as.factor(merged.dt$environment)
  merged.dt$range<-as.factor(merged.dt$range)
    merged.dt$pass<-as.factor(merged.dt$pass)
merged.dt$pedigree<-as.factor(merged.dt$pedigree)
  merged.dt$tester<-as.factor(merged.dt$tester)
            
#subet out the covariate and rerun the spatial analysis
rstructs <-c('ar1v(range):ar1(pass)', 'idv(range):ar1(pass)', 'ar1(range):idv(pass)') #assuming the pass is within the row 
spatial.list<-list()
 
#Test each of the residual structures  
for (r in rstructs){
  if("None" != Covariate){
        mergedFixed<-copy(merged.dt)
#..i) Spatial model as a fixed effect ####
      spatial_mod<-tryCatch(asreml(fixed = pheno_value ~ pedigree + value, random = ~ rep, 
                            residual = as.formula(paste("~",r)),na.action = na.method(x=c("include")),  data = mergedFixed),error=function(x)return(NA))
    }else{ 
#....Model without a covariate
    no.covar.melt<-merge(sub.pheno,r_p,by.x=c("environment","range","pass"), by.y = c("environment","range","pass"),all = T)
      no.covar.melt$range<-as.factor(no.covar.melt$range)
        no.covar.melt$pass<-as.factor(no.covar.melt$pass)
#Model    
    spatial_mod<-tryCatch(asreml(fixed = value ~ pedigree, random = ~ rep, residual =  as.formula(paste("~",r)),
                    na.action = na.method(x="include"),  data = no.covar.melt),error=function(x)return(NA))}
              
#Results of spatial analysis
  if(!is.na(spatial_mod$loglik)){
    spatial_DT<-as.data.table(summary(spatial_mod)$varcomp,keep.rownames="id")
    spatial_DT[id=="value",id:=paste(Covariate)]
    spatial_DT$Covariate<-Covariate
    spatial_DT$BIC<-BIC(spatial_mod)
    spatial_DT$structure<-r
#...ensure the correlations are stable
      cor.components = grep(".cor$", names(spatial_mod$vparameters), value = T)
       max.ar1 = SpatCor #user defined perameter set at 0.60 for default based on Rogers et al. (2021)
        cors_resid<-abs(spatial_mod$vparameters[names(spatial_mod$vparameters) %in% cor.components])
          spatial_DT$Correlation_Instability<-max(cors_resid)
spatial.list[[r]]<-spatial_DT

#Save the values in a list          
  complete.spatial<-rbindlist(spatial.list)
    complete.spatial$tester<-testers[t]
      complete.spatial$pheno<-pheno[i]
  }else{complete.spatial<-NULL}
} #end the first loop of the spatial

#..B) Heritability from spatial analysis ####
if(!is.na(spatial_mod$loglik)){  
  resid.struct<-unique(complete.spatial[, .SD[which.min(BIC)]]$structure)
  RESID<-c("id(units)",resid.struct)
resid_mat<-matrix(nrow = length(RESID), ncol=3, dimnames = list(c(RESID),c("Residual","H2_cullis","Correlation_Instability")))
#...i) Calculate the heritability for each of the spatial models with a covariate ####      
            for (k in 1:length(RESID)){
              if("None" != Covariate){
                  mergedRandom<-copy(merged.dt)
                spatial_mod<-asreml(fixed = pheno_value ~ value, random = ~ pedigree + rep,na.action = na.method(x="include"), residual = as.formula(paste("~",RESID[k])),
                                      data = mergedRandom)}else{
#...ii) Calculate the heritability for each of the spatial models without a covariate ####      
   no.covar.melt<-merge(sub.pheno,r_p,by.x=c("environment","range","pass"), by.y = c("environment","range","pass"),all = T)
#....add the correct case
   no.covar.melt$environment<-as.factor(no.covar.melt$environment)
    no.covar.melt$range<-as.factor(no.covar.melt$range)
      no.covar.melt$pass<-as.factor(no.covar.melt$pass)
    no.covar.melt$pedigree<-as.factor(no.covar.melt$pedigree)
      no.covar.melt$tester<-as.factor(no.covar.melt$tester)
        no.covar.meltCopy<-copy(no.covar.melt)
# Random Spatial Model without a covariate 
  spatial_mod<-asreml(fixed = value~ 1, random = ~ pedigree + rep, residual =  as.formula(paste("~",RESID[k])),
                            na.action = na.method(x="include"),  data = no.covar.meltCopy)}
 #...Save the Spatial analysis terms                     
    pred_h2<-as.data.table(summary(spatial_mod)$varcomp,keep.rownames="id")
      pred_spat<-predict.asreml(spatial_mod,classify = "pedigree", present = "pedigree")
        resid_mat[k,"H2_cullis"]<-(1 - ((pred_spat$avsed**2) / (2 * pred_h2[id=="pedigree"]$component)))* 100
        resid_mat[k,"Residual"]<-RESID[k]
              
#...ensure the correlations are stable
    cor.components = grep(".cor$", names(spatial_mod$vparameters), value = T)
        max.ar1 = SpatCor
    if (length(cor.components > 0)){
                cors_resid<-abs(spatial_mod$vparameters[names(spatial_mod$vparameters) %in% cor.components])
                  rm(spatial_mod)
                resid_mat[k,"Correlation_Instability"]<-max(cors_resid)
      }else{resid_mat[k,"Correlation_Instability"]<-NA}
    }
#Table for heritability after spatial analysis
    h2_resid.DT<-as.data.table(resid_mat)
    h2_resid.DT$H2_cullis<-as.numeric(h2_resid.DT$H2_cullis)
    h2_resid.DT$Correlation_Instability<-as.numeric(h2_resid.DT$Correlation_Instability)
#..add the phenotype, tester, and environment then store the results  
      h2_resid.DT[,`:=`(pheno=paste(pheno[i]), tester=paste(testers[t]), environment=min(locs[j]),Covariate=paste(Covariate))]          
    }
}else{complete.spatial<-NULL} #ends the spatial analysis.... else statement to still generate a list for scenairos where a spatial analysis will not have taken place 

#Save the BIC of the Models
  ModStats<-Wald.Test_DT[,.SD[which.max(BIC)], by=.(Covariate,pheno,tester,environment)][,c('Covariate','pheno','tester','environment','BIC')]
    ModStats[,Correlation_Instability:=NA]  
        ModStats$Model<-'IID'
      if(!is.null(complete.spatial)){
        SpatMods<-complete.spatial[,.SD[which.max(BIC)], by=.(Covariate,pheno,tester,structure)][,environment:=paste(locs[j])][,c('Covariate','pheno','tester','environment','structure','BIC','Correlation_Instability')]
          setnames(SpatMods,c('structure'),c('Model'))
            ModStats<-rbind(ModStats,SpatMods)
    }  
ModelBIC[[paste0(pheno[i],'_',testers[t],'_',locs[j])]]<-ModStats
      
#3) BLUEs by environment using the model with the best fit ####
r.struct<-ModStats[Correlation_Instability < SpatCor | is.na(Correlation_Instability), .SD[which.min(BIC)]]$Model #will get the best residual model based on BIC
    r.struct<-gsub('IID','units',r.struct)
      Covariate<-sig.Her.loc[, .SD[which.max(H2_cullis)]]$Covariate #identify the proper covariate based on maximum heritability 

#..A) Run Final Model with lowest BIC ####
              
# Model for covariate     
if(Covariate=='None'){     
    if(r.struct !='units'){ #For spatial analysis with no covariate 
      ModFinal<-asreml(fixed = value ~ pedigree, random = ~rep,na.action = na.method(x="include"),residual = as.formula(paste("~",r.struct)), data = merged.dt)
      }else{ #For No spatial analysis with no covariate
        ModFinal<-asreml(fixed = value ~ pedigree, random = ~rep,na.action = na.method(x="include"),residual = as.formula(paste("~",r.struct)), data = sub.pheno)
  }
}else{
  if(r.struct !='units'){ #For spatial analysis with a covariate
    ModFinal<-asreml(fixed = pheno_value ~ pedigree + value, random = ~ rep,
                     na.action = na.method(x="include"), residual = as.formula(paste("~",r.struct)),data = merged.dt)
  }else{ #For No spatial analysis with a covariate
    ModFinal<-asreml(fixed = pheno_value~ pedigree + value, random = ~ rep, 
                     residual = as.formula(paste("~",r.struct)),na.action = na.method(x=c("include")), data = w10004_sub_melt[variable==Covariate])
  }
}
  
pred_vals<-predict.asreml(ModFinal, classify = "pedigree",vcov = T)
#..Save the values              
  Complete.Pheno<-as.data.table(pred_vals$pvals)
    Complete.Pheno<-Complete.Pheno[!is.na(predicted.value)]
      Complete.Pheno$S2.weight<-1 / (sqrt(Complete.Pheno$std.error))
#..Variance components                
    VC<-as.data.table(summary(ModFinal)$varcomp,keep.rownames = "id")
     wald_stat<-as.data.table(wald(ModFinal),keep.rownames = "id")
       wald_stat[id=="value",id:=Covariate]

#.. variance covariance 
   VcovList[[paste0(pheno[i],'_',testers[t],'_',locs[j])]]<-pred_vals$vcov
  
#4) Save the Results ####
VC[,`:=`(pheno=paste(pheno[i]), tester=paste(testers[t]), environment=paste(locs[j]))]
  VarCompList[[paste0(pheno[i],'_',testers[t],'_',locs[j])]]<-VC
wald_stat[,`:=`(pheno=paste(pheno[i]), tester=paste(testers[t]), environment=paste(locs[j]))]
  WaldTestList[[paste0(pheno[i],'_',testers[t],'_',locs[j])]]<-wald_stat
  CovList[[paste0(pheno[i],'_',testers[t],'_',locs[j])]]<-Wald.Test_DT  
Complete.Pheno[,`:=`(pheno=paste(pheno[i]), tester=paste(testers[t]), environment=paste(locs[j]))]
  BLUEsList[[paste0(pheno[i],'_',testers[t],'_',locs[j])]]<-Complete.Pheno
  H2CovariateList[[paste0(pheno[i],'_',testers[t],'_',locs[j])]]<-Her.loc[,.SD[which.max(H2_cullis)],by=.(Covariate)][,c('Covariate','environment','tester','pheno','H2_cullis','BIC')]
  
if(!is.null(complete.spatial)){   #in case no environment has range and pass information or the model failed 
    h2_resid.DT[Residual=='id(units)',Residual:='IID']
  H2ResidualList[[paste0(pheno[i],'_',testers[t],'_',locs[j])]]<-h2_resid.DT
}else{
  Her.loc[,Correlation_Instability:=NA]
    Her.loc[,Residual:='IID']
  H2ResidualList[[paste0(pheno[i],'_',testers[t],'_',locs[j])]]<- Her.loc[,.SD[which.min(BIC)]][,c('Residual','H2_cullis','Correlation_Instability','pheno','tester','environment','Covariate')]
}
print(paste0(pheno[i],'_',testers[t],'_',locs[j], ' Complete!'))
    }
  }
print(paste0(locs[j],' is done ',round(j/length(locs)*100,digits = 2),'% Environments complete'))
}
return(list(VarCompList,WaldTestList,CovList,BLUEsList,H2ResidualList,H2CovariateList,VcovList,ModelBIC,skip.list))}

