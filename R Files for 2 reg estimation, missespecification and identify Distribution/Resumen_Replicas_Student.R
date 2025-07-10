####Resume N replications----
#rm(list = ls())

###load simulations for estimation, you can choose estimation with the true distribution or with any other, for instance
### gaussian to get respective column of table 1 and also complete table 10,  this is the example for True distribution error student-t  for M1 model
####If you have only one file con all the results of the simulations use:
load("replicas_Identify_Dist_RelativeBias_Forecasting_1000_2reg_ForStudent.rds")

###If you need to load each file and concatenate the outputs of the simulations use, however in this case that lines are with comment.
###
#load("replicas_Identify_Dist_RelativeBias_Forecasting_1000_2reg_ForStudent.rds"")

#repl_estimation<-repl_student_estimation
#repl_compare_dist_todo<-repl_student_compare_dist
#repl<-repl_student

###Concatenate several files for Estimation and Credible Intervals 

#repl_estimation<-c(repl_estimation,repl_student_estimation)

###Concatenate several files for simulations DIC and WAIC for compare distributions, name of the 

#repl_compare_dist_todo<-c(repl_compare_dist_todo,repl_student_compare_dist)

##Concatenate several files for Forecasting n.ahead with TRUE distribution

#repl<-c(repl,repl_student)





any(sapply(repl_estimation, is.null)) ##Verify null entries
any(sapply(repl, is.null)) ##Verify null entries
any(sapply(repl_compare_dist_todo, is.null)) ##Verify null entries


 
n_rep=1000 ###Set number of simulations


fun = function(x){x/n_rep}



regimes=2
veces_est<-vector("list",regimes)
longitud_est<-vector("list",regimes)
sd_est<-vector("list",regimes)



h.ahead=10
k = 3 ##Dimension of the output vector
v=2 ## Dimension of the exogenous variable vector
ars <- list(p=c(1,2),q=c(1,0),d=c(1,0))


###Create vectors for get counts
veces<-vector("list",h.ahead)
longitud<-vector("list",h.ahead)
sd<-vector("list",h.ahead)




###Counts for estimation ----
for(reg in 1:regimes){
  dime_inter<-dim(Rg[[reg]]$cs)
  dime_phi=list()
  for(orders in 1:ars$p[reg]){
       dime_phi<-append(dime_phi,list(dim(Rg[[reg]]$Location$phi1)))
     }
  dime_sigma<-dim(Rg[[reg]]$Sigma)
  suma_intercept<-matrix(rep(0,prod(dime_inter)),dime_inter[1],dime_inter[2])
  suma_phi=list()
  for(orders in 1:ars$p[reg]){
    suma_phi<-append(suma_phi,list(matrix(rep(0,prod(dime_phi[[orders]])),dime_phi[[orders]][1],dime_phi[[orders]][2])))
  }
  if(reg==1)
  {
    dime_beta=list()
    for(orders in 1:ars$q[reg]){
      dime_beta<-append(dime_beta,list(dim(Rg[[reg]]$Location$beta1)))}
    dime_delta=list()
    for(orders in 1:ars$d[reg]){
      dime_delta<-append(dime_delta,list(dim(Rg[[reg]]$Location$delta1)))}
  
    
    suma_beta=list()
    for(orders in 1:ars$p[reg]){
      suma_beta<-append(suma_beta,list(matrix(rep(0,prod(dime_beta[[orders]])),dime_beta[[orders]][1],dime_beta[[orders]][2])))}
      
      suma_delta=list()
      for(orders in 1:ars$d[reg]){
        suma_delta<-append(suma_delta,list(matrix(rep(0,prod(dime_delta[[orders]])),dime_delta[[orders]][1],dime_delta[[orders]][2])))}
  } 
  
  
  suma_sigma<-matrix(rep(0,prod(dime_sigma)),dime_sigma[1],dime_sigma[2])
  for(rep in 1:n_rep){
    
    ###Autoregressive
    row.names.obj_Reg_autoreg<-rownames(repl_estimation[[rep]]$non_structural[[reg]]$Autoregressive)
    col.names.obj_Reg_autoreg<-colnames(repl_estimation[[rep]]$non_structural[[reg]]$Autoregressive)
    pos.mean<-which(colnames(repl_estimation[[rep]]$non_structural[[reg]]$Autoregressive)==col.names.obj_Reg_autoreg[1])
    pos.low<-which(colnames(repl_estimation[[rep]]$non_structural[[reg]]$Autoregressive)==col.names.obj_Reg_autoreg[4])
    pos.high<-which(colnames(repl_estimation[[rep]]$non_structural[[reg]]$Autoregressive)==col.names.obj_Reg_autoreg[5])
    pos.sd<-which(colnames(repl_estimation[[rep]]$non_structural[[reg]]$Autoregressive)==col.names.obj_Reg_autoreg[2])
    
    ####IC
    ###low
    int.low_Reg<-repl_estimation[[rep]]$non_structural[[reg]]$Autoregressive[,pos.low]
    int.low.intercep_Reg<-int.low_Reg[row.names.obj_Reg_autoreg[1],]
    int.low_Phi_Reg<-t(int.low_Reg[row.names.obj_Reg_autoreg[2:(k+1)],])
    
      
      
    if(ars$q[reg]!=0){
    int.low_Beta_Reg<-t(int.low_Reg[row.names.obj_Reg_autoreg[(k+2):(k+1+v)],])
    }
    
    if(ars$d[reg]!=0){
      int.low_Delta_Reg<-t(int.low_Reg[row.names.obj_Reg_autoreg[(k+1+v+1):(k+1+v+1)],])
    }
    ###high
    int.high_Reg<-repl_estimation[[rep]]$non_structural[[reg]]$Autoregressive[,pos.high]
    int.high.intercep_Reg<-int.high_Reg[row.names.obj_Reg_autoreg[1],]
    int.high_Phi_Reg<-t(int.high_Reg[row.names.obj_Reg_autoreg[2:(k+1)],])
    if(ars$p[reg]>1){
      int.low_Phi2_Reg<-t(int.low_Reg[row.names.obj_Reg_autoreg[(k+2):(2*k+1)],])
      int.high_Phi2_Reg<-t(int.high_Reg[row.names.obj_Reg_autoreg[(k+2):(2*k+1)],])
    }
    
    if(ars$q[reg]!=0){
    int.high_Beta_Reg<-t(int.high_Reg[row.names.obj_Reg_autoreg[(k+2):(k+1+v)],])
    }
    
    if(ars$d[reg]!=0){
      int.high_Delta_Reg<-t(int.high_Reg[row.names.obj_Reg_autoreg[(k+1+v+1):(k+1+v+1)],])
    }
    
    
    ####Phi
    
    
    Phi_TF<-(int.low_Phi_Reg<Rg[[reg]]$Location$phi1) & (Rg[[reg]]$Location$phi1<int.high_Phi_Reg)
    Include_Phi<- matrix(as.integer(as.logical(Phi_TF)),dime_phi[[reg]][1],dime_phi[[reg]][2])
    suma_phi[[1]]<-suma_phi[[1]]+Include_Phi
    
    if(reg==2){
      Phi2_TF<-(int.low_Phi2_Reg<Rg[[reg]]$Location$phi2) & (Rg[[reg]]$Location$phi2<int.high_Phi2_Reg)
      Include_Phi2<- matrix(as.integer(as.logical(Phi2_TF)),dime_phi[[reg]][1],dime_phi[[reg]][2])
      suma_phi[[2]]<-suma_phi[[2]]+Include_Phi2  
    }
    
    
    ####Comparación
    Intercept_TF<-(int.low.intercep_Reg<Rg[[reg]]$cs) & (Rg[[reg]]$cs<int.high.intercep_Reg)
    Include_Intercept<- matrix(as.integer(as.logical(Intercept_TF)),dime_inter[1],dime_inter[2])
    suma_intercept<-suma_intercept+Include_Intercept
    
    
    
    ####Beta
    if(ars$q[reg]!=0){
    Beta_TF<-(int.low_Beta_Reg<Rg[[reg]]$Location$beta1) & (Rg[[reg]]$Location$beta1<int.high_Beta_Reg)
    Include_Beta<- matrix(as.integer(as.logical(Beta_TF)),dime_beta[[reg]][1],dime_beta[[reg]][2])
    suma_beta[[reg]]<-suma_beta[[reg]]+Include_Beta
    }
    
    ####Delta
    if(ars$d[reg]!=0){
      Delta_TF<-(int.low_Delta_Reg<t(Rg[[reg]]$Location$delta1)) & (t(Rg[[reg]]$Location$delta1)<int.high_Delta_Reg)
      Include_Delta<- matrix(as.integer(as.logical(Delta_TF)),dime_delta[[reg]][1],dime_delta[[reg]][2])
      suma_delta[[reg]]<-suma_delta[[reg]]+Include_Delta
    }
    
    ###Sigma
    row.names.obj_Reg_scale<-rownames(repl_estimation[[rep]]$non_structural[[reg]]$Scale)
    col.names.obj_Reg_scale<-colnames(repl_estimation[[rep]]$non_structural[[reg]]$Scale)
    pos.mean_scale<-which(colnames(repl_estimation[[rep]]$non_structural[[reg]]$Scale)==col.names.obj_Reg_scale[1])
    pos.low_scale<-which(colnames(repl_estimation[[rep]]$non_structural[[reg]]$Scale)==col.names.obj_Reg_scale[2])
    pos.high_scale<-which(colnames(repl_estimation[[rep]]$non_structural[[reg]]$Scale)==col.names.obj_Reg_scale[3])
    
    ####IC
    ###low
    int.low_Reg_scale<-repl_estimation[[rep]]$non_structural[[reg]]$Scale[,pos.low_scale]
    
    int.low_Sigma_Reg<-t(int.low_Reg_scale)
    
    ###high
    int.high_Reg_scale<-repl_estimation[[rep]]$non_structural[[reg]]$Scale[,pos.high_scale]
    
    int.high_Sigma_Reg<-t(int.high_Reg_scale)
    
    
    ####Comparación
    
    Sigma_TF<-(int.low_Sigma_Reg<(Rg[[reg]]$Sigma)) & ((Rg[[reg]]$Sigma)<int.high_Sigma_Reg)
    Include_Sigma<- matrix(as.integer(as.logical(Sigma_TF)),dime_sigma[1],dime_sigma[2])
    suma_sigma<-suma_sigma+Include_Sigma
    
    print(rep)
  }
  if(reg==1){
  veces_est[[reg]]<-list(Intercept=suma_intercept,Phi=suma_phi[[1]],Beta=suma_beta[[1]],Delta=suma_delta[[1]],Sigma=suma_sigma)
  }
  if(reg==2)
  {veces_est[[reg]]<-list(Intercept=suma_intercept,Phi=suma_phi[[1]],Phi2=suma_phi[[2]],Sigma=suma_sigma)}
}

#para.extra=FALSE
para.extra=TRUE
if(isTRUE(para.extra)){

extra_1<-extra
suma_extra<-rep(0,length(extra_1))
names.extra<-colnames(repl_estimation[[1]]$extra)
pos.low.extra<-which(colnames(repl_estimation[[1]]$extra)==names.extra[4])
pos.high.extra<-which(colnames(repl_estimation[[1]]$extra)==names.extra[5])
sesgo_extra<-rep(0,length(extra_1))
}


suma_thresholds<-matrix(rep(0,regimes-1),regimes-1,1)
dife_thresholds<-matrix(rep(0,regimes-1),regimes-1,1)
r_true<-as.matrix(umbrales)
delay<-rep(0,n_rep)

for(rep in 1:n_rep){
  dife_thresholds=dife_thresholds+(repl_estimation[[rep]]$thresholds[,1]-r_true)
  delay[rep]<-repl_estimation[[rep]]$d
  
  if(isTRUE(para.extra)){
    for(l in 1:length(extra_1)){
  	sesgo_extra=sesgo_extra+(repl_estimation[[rep]]$extra[,l]-extra_1)
  
  if((repl_estimation[[rep]]$extra[l,pos.low.extra]<extra_1[l]) & ((repl_estimation[[rep]]$extra[l,pos.high.extra]>extra_1[l])) )
  {suma_extra[l]=suma_extra[l]+1}
    }
  }
  
  Threshold_TF<-r_true>repl_estimation[[rep]]$thresholds[,2] & r_true<repl_estimation[[rep]]$thresholds[,3]
  Include_Threshold<-matrix(as.integer(as.logical(Threshold_TF)),regimes-1,1)
  suma_thresholds=suma_thresholds+Include_Threshold
}

sesgo_thresholds<-dife_thresholds/n_rep
#sesgo_thresholds

divided<-function(x){x/10}


prop_Reg1<-lapply(veces_est[[1]],fun)
prop_Reg2<-lapply(veces_est[[2]],fun)
#prop_Reg3<-lapply(veces_est[[3]],fun)

###Percentage of time credible intervals captures real individual parameters table 1 respective column 
lapply(veces_est[[1]],divided)
lapply(veces_est[[2]],divided)

(table(delay)/n_rep)*100 ##Porcentage of times  delay is selected correctly
#suma_extra/n_rep
prop.thresholds<-suma_thresholds/n_rep
prop.delay<-table(delay)/n_rep

prop.thresholds*100
prop.delay*100


if(para.extra){
prop.extra=suma_extra/n_rep
sesgo_extra_def<-sesgo_extra/n_rep
prop.extra*100
sesgo_relativo_extra<-(sesgo_extra_def/abs(extra))*100
sesgo_relativo_extra
}



##Relative Bias: extra and thresholds get Table 2

sesgo_relativo_thresholds<-(sesgo_thresholds/abs(umbrales))*100


sesgo_relativo_thresholds





#####forecasting----

#repl<-repl_laplace
#repl<-c(repl,repl_student)
for(pasos in 1:h.ahead){
  suma_pasos=rep(0,k)
  suma_longitud=rep(0,k)
  #suma_sd=rep(0,k)
  for(rep in 1:n_rep){
    int.low.pos=which(names(repl[[rep]]$forecasting[pasos,])=="HDI_Low")
    int.high.pos=which(names(repl[[rep]]$forecasting[pasos,])=="HDI_high")
    #sd.pos=which(names(repl[[rep]]$forecasting[pasos,])=="SD")
    for(entrada in 1:k){
      if((repl[[rep]]$True_Values[pasos,entrada]>repl[[rep]]$forecasting[pasos,int.low.pos[entrada]])&(repl[[rep]]$True_Values[pasos,entrada]<repl[[rep]]$forecasting[pasos,int.high.pos[entrada]])){suma_pasos[entrada]=suma_pasos[entrada]+1}
      suma_longitud[entrada]<-suma_longitud[entrada]+(repl[[rep]]$forecasting[pasos,int.high.pos[entrada]]-repl[[rep]]$forecasting[pasos,int.low.pos[entrada]])
      #suma_sd[entrada]<-suma_sd[entrada]+repl[[rep]]$forecasting[pasos,sd.pos[entrada]]
    }
    #print(suma_longitud_gauss)
    #print(suma_sd_gauss)
  }
  veces[[pasos]]=suma_pasos
  longitud[[pasos]]=suma_longitud
  #sd[[pasos]]=suma_sd
}

mult<-function(x){x*100}
prop_pred_int<-lapply(veces,fun)
longitud_med_pred_int_student<-lapply(longitud,fun)

###column of table 5
print(lapply(prop_pred_int,mult),digits=3)



prop_pred_int_student<-prop_pred_int
longitud_med_pred_int_student<-longitud_med_pred_int_student



#load("replicas_gsl.rds")
#save(prop_student_estimation,prop.delay,prop.thresholds,sesgo_thresholds,prop.extra,sesgo_extra_def,file="resume_estimation_student10002reg.rds")
#save(prop_pred_int_student,longitud_med_pred_int_student,file="resume_prediction_student10002reg.rds")
     
#######Relative Bias for missespecification distribution----
####True distribution the:
#repl_estimation_oth_dist<-repl_estimation
repl_estimation_oth_dist<-repl_estimation_oth_dist_laplace

any(sapply(repl_estimation_oth_dist, is.null))

########Optional if you need to concatenate some files
repl_estimation_oth_dist_Gaussian<-repl_gaussian_estimation
repl_estimation_oth_dist_slash<-repl_slash_estimation
repl_estimation_oth_dist_contaminated<-repl_contaminated_estimation
repl_estimation_oth_dist_hyperbolic<-repl_hyperbolic_estimation
repl_estimation_oth_dist_laplace<-repl_laplace_estimation


#repl_estimation_oth_dist_slash<-c(repl_estimation_oth_dist_slash,repl_slash_estimation)
#repl_estimation_oth_dist_contaminated<-c(repl_estimation_oth_dist_contaminated,repl_contaminated_estimation)
#repl_estimation_oth_dist_hyperbolic<-c(repl_estimation_oth_dist_hyperbolic,repl_hyperbolic_estimation)
#repl_estimation_oth_dist_laplace<-c(repl_estimation_oth_dist_laplace,repl_laplace_estimation)
#repl_estimation_oth_dist_Gaussian<-c(repl_estimation_oth_dist,repl_gaussian_estimation)





n_rep=1000
fun = function(x){x/n_rep}
veces<-vector("list",h.ahead)
longitud<-vector("list",h.ahead)
sd<-vector("list",h.ahead)

regimes=2
suma_est<-vector("list",regimes)
longitud_est<-vector("list",regimes)
sd_est<-vector("list",regimes)


###Estimación
for(reg in 1:regimes){
  dime_inter<-dim(Rg[[reg]]$cs)
  dime_phi=list()
  for(orders in 1:ars$p[reg]){
    dime_phi<-append(dime_phi,list(dim(Rg[[reg]]$Location$phi1)))
  }
  dime_sigma<-dim(Rg[[reg]]$Sigma)
  suma_intercept<-matrix(rep(0,prod(dime_inter)),dime_inter[1],dime_inter[2])
  suma_phi=list()
  for(orders in 1:ars$p[reg]){
    suma_phi<-append(suma_phi,list(matrix(rep(0,prod(dime_phi[[orders]])),dime_phi[[orders]][1],dime_phi[[orders]][2])))
  }
  if(reg==1)
  {
    dime_beta=list()
    for(orders in 1:ars$q[reg]){
      dime_beta<-append(dime_beta,list(dim(Rg[[reg]]$Location$beta1)))}
    dime_delta=list()
    for(orders in 1:ars$d[reg]){
      dime_delta<-append(dime_delta,list(dim(Rg[[reg]]$Location$delta1)))}
    
    
    suma_beta=list()
    for(orders in 1:ars$p[reg]){
      suma_beta<-append(suma_beta,list(matrix(rep(0,prod(dime_beta[[orders]])),dime_beta[[orders]][1],dime_beta[[orders]][2])))}
    
    suma_delta=list()
    for(orders in 1:ars$d[reg]){
      suma_delta<-append(suma_delta,list(matrix(rep(0,prod(dime_delta[[orders]])),dime_delta[[orders]][1],dime_delta[[orders]][2])))}
  } 
  
  
  suma_sigma<-matrix(rep(0,prod(dime_sigma)),dime_sigma[1],dime_sigma[2])
  for(rep in 1:n_rep){
    
    ###Autoregressive
    row.names.obj_Reg_autoreg<-rownames(repl_estimation_oth_dist[[rep]]$non_structural[[reg]]$Autoregressive)
    col.names.obj_Reg_autoreg<-colnames(repl_estimation_oth_dist[[rep]]$non_structural[[reg]]$Autoregressive)
    pos.mean<-which(colnames(repl_estimation_oth_dist[[rep]]$non_structural[[reg]]$Autoregressive)==col.names.obj_Reg_autoreg[1])
    
    ####Estimation
    
    est.Reg<-repl_estimation_oth_dist[[rep]]$non_structural[[reg]]$Autoregressive[,pos.mean]
    est.intercep_Reg<-est.Reg[row.names.obj_Reg_autoreg[1],]
    est_Phi_Reg<-t(est.Reg[row.names.obj_Reg_autoreg[2:(k+1)],])
    if(ars$p[reg]>1){
      est_Phi2_Reg<-t(est.Reg[row.names.obj_Reg_autoreg[(k+2):(2*k+1)],])
    }
    
    
    
    if(ars$q[reg]!=0){
      est_Beta_Reg<-t(est.Reg[row.names.obj_Reg_autoreg[(k+2):(k+1+v)],])
    }
    
    if(ars$d[reg]!=0){
      est_Delta_Reg<-t(est.Reg[row.names.obj_Reg_autoreg[(k+1+v+1):(k+1+v+1)],])
    }
    
    ####Phi
    
    
    
    suma_phi[[1]]<-suma_phi[[1]]+(est_Phi_Reg-Rg[[reg]]$Location$phi1)
    
    if(reg==2){
      suma_phi[[2]]<-suma_phi[[2]]+(est_Phi2_Reg-Rg[[reg]]$Location$phi2)  
    }
    
    
    ####Intercepto
    
    suma_intercept<-suma_intercept+(est.intercep_Reg-Rg[[reg]]$cs)
    
    
    
    ####Beta
    if(ars$q[reg]!=0){
      suma_beta[[reg]]<-suma_beta[[reg]]+(est_Beta_Reg-Rg[[reg]]$Location$beta1)
    }
    
    ####Delta
    if(ars$d[reg]!=0){
      suma_delta[[reg]]<-suma_delta[[reg]]+(t(est_Delta_Reg)-Rg[[reg]]$Location$delta1)
    }
    
    ###Sigma
    row.names.obj_Reg_scale<-rownames(repl_estimation_oth_dist[[rep]]$non_structural[[reg]]$Scale)
    col.names.obj_Reg_scale<-colnames(repl_estimation_oth_dist[[rep]]$non_structural[[reg]]$Scale)
    pos.mean_scale<-which(colnames(repl_estimation_oth_dist[[rep]]$non_structural[[reg]]$Scale)==col.names.obj_Reg_scale[1])
    
    ####IC
    ###low
    est_Reg_scale<-repl_estimation_oth_dist[[rep]]$non_structural[[reg]]$Scale[,pos.mean_scale]
    
    
    ####Comparación
    
    
    suma_sigma<-suma_sigma+(est_Reg_scale-Rg[[reg]]$Sigma)
    
    #print(rep)
  }
  if(reg==1){
    suma_est[[reg]]<-list(Intercept=suma_intercept,Phi1=suma_phi[[1]],Beta=suma_beta[[1]],Delta=suma_delta[[1]],Sigma=suma_sigma)
  }
  if(reg==2)
  {suma_est[[reg]]<-list(Intercept=suma_intercept,Phi1=suma_phi[[1]],Phi2=suma_phi[[2]],Sigma=suma_sigma)}
}

para.extra=FALSE

if(isTRUE(para.extra)){
  
  extra_1<-extra
  suma_extra<-rep(0,length(repl_estimation_oth_dist[[1]]$extra[,1]))
  names.extra<-colnames(repl_estimation_oth_dist[[1]]$extra)
  pos.mean.extra<-which(colnames(repl_estimation_oth_dist[[1]]$extra)==names.extra[1])
  sesgo_extra<-rep(0,length(repl_estimation_oth_dist[[1]]$extra[,1]))
  suma_est_extra<-rep(0,length(repl_estimation_oth_dist[[1]]$extra[,1]))
}


suma_thresholds<-matrix(rep(0,regimes-1),regimes-1,1)
dife_thresholds<-matrix(rep(0,regimes-1),regimes-1,1)
r_true<-as.matrix(umbrales)
delay<-rep(0,n_rep)

for(rep in 1:n_rep){
  dife_thresholds=dife_thresholds+(repl_estimation_oth_dist[[rep]]$thresholds[,1]-r_true)
  delay[rep]<-repl_estimation_oth_dist[[rep]]$d
  
  if(isTRUE(para.extra)){
    for(l in 1:length(extra_1)){
      sesgo_extra=sesgo_extra+(repl_estimation_oth_dist[[rep]]$extra[,l]-extra_1)
      suma_est_extra=suma_est_extra+repl_estimation_oth_dist[[rep]]$extra[,l]
      if((repl_estimation_oth_dist[[rep]]$extra[l,pos.low.extra]<extra_1[l]) & ((repl_estimation_oth_dist[[rep]]$extra[l,pos.high.extra]>extra_1[l])) )
      {suma_extra[l]=suma_extra[l]+1}
    }
  }
  
  Threshold_TF<-r_true>repl_estimation_oth_dist[[rep]]$thresholds[,2] & r_true<repl_estimation_oth_dist[[rep]]$thresholds[,3]
  Include_Threshold<-matrix(as.integer(as.logical(Threshold_TF)),regimes-1,1)
  suma_thresholds=suma_thresholds+Include_Threshold
}

sesgo_thresholds<-dife_thresholds/n_rep


promedio_est_extra<-suma_est_extra/n_rep


sesgo_Reg1<-lapply(suma_est[[1]],fun)
sesgo_Reg2<-lapply(suma_est[[2]],fun)
###Regimen 1
print(sesgo_Reg1$Intercept/Rg[[1]]$cs*100,digits=2)
print(sesgo_Reg1$Phi1/Rg[[1]]$Location$phi1*100,digits=2)
print(sesgo_Reg1$Beta/Rg[[1]]$Location$beta1*100,digits=2)
print(sesgo_Reg1$Delta/Rg[[1]]$Location$delta1*100,digits=2)
print(sesgo_Reg1$Sigma/Rg[[1]]$Sigma*100,digits=4)

###Regimen 2
print(sesgo_Reg2$Intercept/Rg[[2]]$cs*100,digits=2)
print(sesgo_Reg2$Phi1/Rg[[2]]$Location$phi1*100,digits=2)
print(sesgo_Reg2$Phi2/Rg[[2]]$Location$phi2*100,digits=2)
print(sesgo_Reg2$Sigma/Rg[[2]]$Sigma*100,digits=4)


(table(delay)/n_rep)*100 ##Porcentaje de veces

if(para.extra){
  prop.extra=suma_extra/n_rep
  sesgo_extra_def<-sesgo_extra/n_rep
  sesgo_relativo_extra<-(sesgo_extra_def/extra)*100
  sesgo_relativo_extra
}

###Sesgos relativos de los parámetros: extra and thresholds

sesgo_relativo_thresholds<-(sesgo_thresholds/umbrales)*100


sesgo_relativo_thresholds


#### Mejor Modelo ----
#repl_compare_dist_todo<-repl_student_compare_dist
#repl_compare_dist_todo<-c(repl_compare_dist_todo,repl_student_compare_dist)

min_DIC_dist=0
min_WAIC_dist=0
for(j in 1:n_rep)
{if(repl_compare_dist_todo[[j]]$Min_DIC_dist[1]=="fit1_2reg_Student") min_DIC_dist=min_DIC_dist+1
if(repl_compare_dist_todo[[j]]$Min_WAIC_dist[1]=="fit1_2reg_Student") min_WAIC_dist=min_WAIC_dist+1
}
min_DIC_dist/n_rep
min_WAIC_dist/n_rep

###Segundo mejor modelo 2reg
#repl_student_compare_dist[[j]]$WAIC_models_dist[,order(repl_student_compare_dist[[1]]$WAIC_models_dist,decreasing=FALSE),drop=FALSE]

min2do_DIC_dist=0
min2do_WAIC_dist=0
for(j in 1:n_rep)
{if(colnames(repl_compare_dist_todo[[j]]$DIC_models_dist[,order(repl_compare_dist_todo[[j]]$DIC_models_dist,decreasing=FALSE),drop=FALSE])[2]=="fit1_2reg_Student") min2do_DIC_dist=min2do_DIC_dist+1
if(colnames(repl_compare_dist_todo[[j]]$WAIC_models_dist[,order(repl_compare_dist_todo[[j]]$WAIC_models_dist,decreasing=FALSE),drop=FALSE])[2]=="fit1_2reg_Student") min2do_WAIC_dist=min2do_WAIC_dist+1
}
min2do_DIC_dist/n_rep
min2do_WAIC_dist/n_rep


save(repl_estimation,repl_compare_dist_todo,repl,repl_estimation_oth_dist_Gaussian,repl_estimation_oth_dist_slash,repl_estimation_oth_dist_contaminated,repl_estimation_oth_dist_hyperbolic,repl_estimation_oth_dist_laplace,umbrales,extra,Rg,ars,h.ahead,file="replicas_Identify_Dist_RelativeBias_Forecasting_1000_2reg_ForStudent.rds")



