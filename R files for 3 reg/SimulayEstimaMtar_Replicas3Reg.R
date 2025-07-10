####Program for simulating MTAR process and estimating its parameters
###This only includes M2 structure model, that is bivariate a MTAR with 3 regimes

# Required packages====

library(GIGrvg)
library(Formula)
library(Rfast)
library(expm)
library(mtar)
####Definición de parámetros ====
n_rep=2
#repl_student<-vector("list", n_rep)
#repl_gaussian<-vector("list", n_rep)
#repl_laplace<-vector("list", n_rep)

#repl_student_estimation<-vector("list", n_rep)
#repl_gaussian_estimation<-vector("list", n_rep)
#repl_laplace_estimation<-vector("list", n_rep)


#repl_hyperbolic<-vector("list", n_rep)
repl_slash<-vector("list", n_rep)
#repl_contaminated<-vector("list", n_rep)

#repl_hyperbolic_estimation<-vector("list", n_rep)
repl_slash_estimation<-vector("list", n_rep)
#repl_contaminated_estimation<-vector("list", n_rep)
calen=100
long=300
h.ahead=10
Tlen = long+calen+h.ahead




##Parámetros del Proceso Yt ====
k = 2
#v=2
v=0
### R1 regimen ====

library(ltsa)
simulate_arma <- function(ar_params = NULL, ma_params = NULL, constant = 0, n = 100, sigma2 = 1){ 
  
  # Load necessary library
  library(stats)
  
  # Validate input
  if (is.null(ar_params) && is.null(ma_params)) {
    stop("Please specify either AR or MA parameters.")
  }
  
  # Generate innovation series
  innovations <- rnorm(n, mean = 0, sd = sqrt(sigma2))
  
  # Initialize series
  arma_series <- numeric(n)
  
  # Handle initial values for AR components
  if (!is.null(ar_params)) {
    arma_series[1:max(length(ar_params), 1)] <- innovations[1:max(length(ar_params), 1)]
  }
  
  # Simulate the ARMA process
  for (t in (max(length(ar_params), 1) + 1):n) {
    ar_part <- ifelse(is.null(ar_params), 0, sum(ar_params * arma_series[(t - 1):(t - length(ar_params))]))
    ma_part <- ifelse(is.null(ma_params), 0, sum(ma_params * innovations[(t - 1):(t - length(ma_params))]))
    arma_series[t] <- constant + ar_part + ma_part + innovations[t]
  }
  return(arma_series)
}




dist <- "Slash"
#extra <- c(0.012,0.08)
#extra<-c(0.1,0.15)
extra=4
delay <- 1
Intercept <- TRUE
ars <- list(p=c(1,1,1),q=c(0,0,0),d=c(0,0,0))
inic<-calen
phi_ar=0.6
sigma2_inn_ar<-1
constant<-1


Location_R1 = list(phi1 = matrix(c(0.8,0,-0.2,0.5),k,k,byrow = TRUE))

Sigma_R1 = matrix(c(1,0,0,4),k,k,byrow = TRUE)

cs_1=matrix(c(2,1),nrow=k)

R1 = list(orders = list(p = 1,q = 0,d = 0),Location = Location_R1,Sigma = Sigma_R1,cs=cs_1)

###Regimen 2
Location_R2 = list(phi1 = matrix(c(0.3,0,0,-0.6),k,k,byrow = TRUE))

Sigma_R2 = matrix(c(1,0,0,1),k,k,byrow = TRUE)

cs_2=matrix(c(0.4,-2),nrow=k)

R2 = list(orders = list(p = 1,q = 0,d = 0),
          Location = Location_R2 ,Sigma = Sigma_R2,cs=cs_2)

###R3 regimen

Location_R3 = list(phi1 = matrix(c(0.6,0,-0.2,0.8),k,k,byrow = TRUE))



Sigma_R3 = matrix(c(2,0,0,1),k,k,byrow = TRUE)

cs_3=matrix(c(-3,0),nrow=k)




R3 = list(orders = list(p = 1,q = 0,d = 0),
          Location_ = Location_R3  ,Sigma = Sigma_R3,cs=cs_3)




## crea lista de objeto tipo Regime
Rg = list(R1 = R1,R2 = R2,R3=R3)


####Mean and Variace ARMA
mean_arma=constant/(1-sum(phi_ar))

var_arma<-ltsa::tacvfARMA(phi=phi_ar,sigma2=sigma2_inn_ar)[1]

params <- list()
umbrales<-qnorm(p=c(0.33,0.66),mean =mean_arma,sd = sqrt(var_arma) )

for(i in 1:length(ars$p)){
  np <- Intercept + ars$p[i]*k + ars$q[i]*v + ars$d[i]
  params[[i]] <- list()
  params[[i]]$location <-rbind(t(Rg[[i]][[4]]),matrix(unlist(Rg[[i]][[2]]),ncol=k,byrow=TRUE))
  #params[[i]]$location <- matrix(c(rbeta(np*k,shape1=4,shape2=16)),np,k)
  # params[[i]]$scale <- diag(rgamma(k,shape=1,scale=1))
  params[[i]]$scale <- Rg[[i]][[3]]
  params[[i]]$scale2 <- chol(params[[i]]$scale)
}
params


## crea lista de objeto tipo Regime

#r = qnorm(p = c(0.33,0.66),mean = 7.1428571,sd=sqrt(2))

#####Alternativa Parámetros#####
#params <- list()

#for(i in 1:length(ars$p)){
#  np <- Intercept + ars$p[i]*k + ars$q[i]*v + ars$d[i]
#  params[[i]] <- list()
#  set.seed(32*i)
#  params[[i]]$location <- round(matrix(c(ifelse(runif(np*k)<=0.5,1,-1)*rbeta(np*k,shape1=4,shape2=16)),np,k),2)
#  set.seed(432*i)
#  params[[i]]$scale <- round(diag(rgamma(k,shape=1,scale=1)),2)
#  params[[i]]$scale2 <- chol(params[[i]]$scale)
#}



start.time <- Sys.time()
bar <- txtProgressBar(min=0, max=n_rep, initial=0, char="=", style=3)
for(replicas in 1:n_rep){
  ## 3 dimensiones ====
  ##2 regimenes con exógenas ====
  ## Obtener el proceso Ut=(Xt,Zt)
  
  #Ut = mtarsim(N = Tlen,Rg = R_ut,seed = NULL)
  #Zt = Ut$Sim$Yt[,3]
  #Xt= Ut$Sim$Yt[,1:2] # Procesos Bidimensional v=2
  
  Zt<-tseries  <- simulate_arma(ar_params = phi_ar, constant = constant, sigma2 = sigma2_inn_ar,n = Tlen+max(ars$p,ars$q,ars$d,delay))
  
  myseries <- matrix(0,Tlen+max(ars$p,ars$q,ars$d,delay),k)
  
  regimen <- cut(Zt[(max(ars$p,ars$q,ars$d,delay)+1-delay):(length(Zt)-delay)],breaks=c(-Inf,umbrales,Inf),labels=1:length(ars$p))
  myseries[1:max(ars$p,ars$q,ars$d),] <- rnorm(max(ars$p,ars$q,ars$d)*k)
  
  for(i in 1:Tlen){
    current <- max(ars$p,ars$q,ars$d,delay) + i
    regimeni <- regimen[i]
    if(Intercept) X <- 1 else X <-  vector()
    for(j in 1:ars$p[regimeni]) X <- c(X,myseries[current-j,])
    if(ars$q[regimeni] > 0) for(j in 1:ars$q[regimeni]) X <- c(X,Xt[current-j,])
    if(ars$d[regimeni] > 0) for(j in 1:ars$d[regimeni]) X <- c(X,Zt[current-j,])
    Theta <- params[[regimeni]]$location   
    mu <- apply(matrix(X,nrow(Theta),ncol(Theta))*Theta,2,sum)
    u <- 1
    if(dist=="Student-t")  u <- 1/rgamma(1,shape=extra/2,rate=extra/2)
    if(dist=="Slash")  u <- 1/rbeta(1,shape1=extra/2,shape2=1)
    if(dist=="Contaminated normal")  if(runif(1)<=extra[1]) u <- 1/extra[2]
    if(dist=="Laplace")  u <- rexp(1,rate=1/8)
    if(dist=="Hyperbolic") u <- rgig(n=1,lambda=1,chi=1,psi=extra^2)
    myseries[current,] <- t(params[[regimeni]]$scale2)%*%matrix(rnorm(k,mean=0,sd=sqrt(u)),k,1) + matrix(mu,k,1)
  }
  
  datos <- data.frame(myseries[(max(ars$p,ars$q,ars$d,delay)+inic+1):length(Zt),])
  colnames(datos) <- paste("X",1:k,sep="")
  datos <- data.frame(datos,threshold=Zt[(max(ars$p,ars$q,ars$d,delay)+inic+1):length(Zt)])
  if(max(ars$q)>0){
    colnames(Xt) <- paste("X",(k+1):(k+v),sep="")
    datos <- data.frame(datos,Xt[(max(ars$p,ars$q,ars$d,delay)+inic+1):length(Zt),])
  }
  #plot(as.ts(datos[,1:k]))
  
  ### Creación de las fechas para datos simulados====
  Fechas=seq(as.Date("2000/1/1"), by = "day", length.out = (Tlen-calen))
  
  
  #### 3 reg 2 dimensiones datos ====
  
  datos1=data.frame(datos,Fecha=Fechas)
  
  
  
  ##  Procedimiento de Estimación Bayesiano  3 reg====
  fecha_final<-Fechas[long]
  
  ##########Student-t or Hyperbolic
  
  fit1_Bayes <- mtar(~X1+X2|threshold, data=datos1, row.names=Fecha,subset={Fecha<=fecha_final} ,ars=ars, dist=dist, n.burnin=500, n.sim=1500, n.thin=1)
  summary_Bayes<-summarymtar_simulation(fit1_Bayes,Get.results = TRUE,Print.Results = FALSE)
  #summarymtar(fit1_Bayes)
  #plot(as.ts(t(fit1_Bayes$chains$extra)))
  
  nano <- forecasting(fit1_Bayes,subset(datos1,Fecha > fecha_final),row.names=Fecha)
  forecasting_salida<-nano$summary
  repl_slash[[replicas]]=list(forecasting=forecasting_salida,True_Values=subset(datos1,Fecha > fecha_final))
  repl_slash_estimation[[replicas]]<-summary_Bayes
  
  
  setTxtProgressBar(bar,replicas)
  #print(replicas)
}

end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken
####Guardar simulaciones====
save(repl_slash_estimation,repl_slash,extra,umbrales,Rg,file="replicas_slash_1000_reg.rds")
