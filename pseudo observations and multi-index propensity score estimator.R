
rm(list = ls())

##load package
library(pseudo)
library(AMORE)

##################################################
##### calculate the true SPCE at time point 5 ####
##################################################
seed <- 2022
set.seed(seed)
n.sample <- 1000000 
rho=sqrt(0.2)
nrho=sqrt(1-rho^2)
Z0=rnorm(n.sample)
Z1=pmin(4,pmax(-4,rnorm(n.sample)*nrho+rho*Z0))
Z2=pmin(4,pmax(-4,rnorm(n.sample)*nrho+rho*Z0))
Z3=pmin(4,pmax(-4,rnorm(n.sample)))
T1=rexp(n.sample,exp(-3-1+0.5*(Z1-Z2)))
T0=rexp(n.sample,exp(-3+0+0.5*(Z1-Z2))) 
est_time <- 5 
Y1_r <- sum(T1 >= est_time)/n.sample
Y0_r <- sum(T0 >= est_time)/n.sample
Real_SPCE <- Y1_r - Y0_r##the true SPCE is 0.144847

##################################
## generate the simulation data ##
##################################
set.seed(seed)
n.sample <- 200 
rho=sqrt(0.2)
nrho=sqrt(1-rho^2)
Z0=rnorm(n.sample)
Z1=pmin(4,pmax(-4,rnorm(n.sample)*nrho+rho*Z0))
Z2=pmin(4,pmax(-4,rnorm(n.sample)*nrho+rho*Z0))
Z3=pmin(4,pmax(-4,rnorm(n.sample)))
A <- rbinom( n.sample, 1, 1/(1+exp(-(Z1+Z2+Z3)/3)) )
TT=rexp(n.sample,exp(-3-A+0.5*(Z1-Z2)))
C=rexp(n.sample,exp(-3.5))
ttime=pmin(TT,C)
event=1*(TT<C)
data <- as.data.frame(cbind(Z1, Z2, Z3, A, ttime, event))


######################################################################
########### the function for calculating the multi-indexes ###########
######################################################################
index.calc <- function(data.pseu,
                       ps1.cov=ps1.cov,ps2.cov=ps2.cov,
                       or1.cov=or1.cov,or2.cov=or2.cov) {
  ##calculate the index from the propensity score model 1
  Si.ps1 <- NA
  if (!is.null(ps1.cov)) {
    form <- as.formula(paste("A~", paste(ps1.cov,collapse="+"), sep = ""))
    fit.ps <- glm(form,data = data.pseu,family = binomial())
    Si.ps1 <- as.matrix(data.pseu[,which(names(data.pseu) %in% ps1.cov)]) %*% fit.ps$coefficients[-1]
    Si.ps1 <- as.numeric(Si.ps1)
  }
  ##calculate the index from the propensity score model 2
  Si.ps2 <- NA
  if (!is.null(ps2.cov)) {
    form <- as.formula(paste("A~", paste(ps2.cov,collapse="+"), sep = ""))
    fit.ps <- glm(form,data = data.pseu,family = binomial())
    Si.ps2 <- as.matrix(data.pseu[,which(names(data.pseu) %in% ps2.cov)]) %*% fit.ps$coefficients[-1]
    Si.ps2 <- as.numeric(Si.ps2)
  }
  
  ##calculate the index from the outcome regression model 1
  Si.rp1 <- NA
  if (!is.null(or1.cov)) {
    form <- as.formula(paste("pseudo~", paste(or1.cov,collapse="+"), "+A",sep = ""))
    beta1=geese(form,scale.fix=TRUE,data=data.pseu,family=gaussian, id=id, jack=F, mean.link="cloglog",corstr="independence")$beta
    Si.rp1 <- as.matrix(data.pseu[,which(names(data.pseu) %in% or1.cov)]) %*% beta1[-c(1,length(beta1))]
    Si.rp1=as.numeric(Si.rp1)
  }
  ##calculate the index from the outcome regression model 2
  Si.rp2 <- NA
  if (!is.null(or2.cov)) {
    form <- as.formula(paste("pseudo~", paste(or2.cov,collapse="+"), "+A",sep = ""))
    beta1=geese(form,scale.fix=TRUE,data=data.pseu,family=gaussian, id=id, jack=F, mean.link="cloglog",corstr="independence")$beta
    Si.rp2 <- as.matrix(data.pseu[,which(names(data.pseu) %in% or2.cov)]) %*% beta1[-c(1,length(beta1))]
    Si.rp2=as.numeric(Si.rp2)
  }
  
  multiindex <- cbind(Si.ps1, Si.ps2,Si.rp1, Si.rp2)
  multiindex <- as.data.frame(multiindex[,!is.na(multiindex[1,])])
  return(multiindex)
}

############################################################################
## function for calculating the SPCE based on peudo-observations and MiPS ##
############################################################################
SPCE.MiPS <- function(data.sim,
                      ps1.cov, ps2.cov, or1.cov, or2.cov,
                      est_time,##the prespecified time point
                      hidden.neurons, learning.rate.global, momentum.global##the hyperparameters
                      ) {
  
  pseudo=pseudosurv(time=data.sim$ttime, event=data.sim$event,tmax=est_time)$pseudo
  id <- c(1:n.sample)
  data <- cbind(data.sim,pseudo,id)
  
  multiindex <- index.calc(data.pseu = data,
                           ps1.cov=ps1.cov,ps2.cov=ps2.cov,
                           or1.cov=or1.cov,or2.cov=or2.cov)
  
  Ai <- as.numeric(data$A)
  n.neurons=c(ncol(multiindex), hidden.neurons, 1)
  net <- newff(n.neurons, learning.rate.global, momentum.global,
               error.criterium = "LMS", Stao = NA, hidden.layer = "tansig", 
               output.layer = "purelin", method = "ADAPTgdwm")
  fit <- train(net, multiindex, Ai, error.criterium="LMS", report=FALSE,
               show.step = 100, n.shows = 5 )
  MiPS <- sim(fit$net,multiindex)
  
  ##calculate the estimated SPCE
  pseudo=data$pseudo
  weight <- data$A/MiPS + (1-data$A)/(1-MiPS)
  surv_1 <- sum(weight[data$A == 1] * pseudo[data$A == 1])/sum(weight[data$A == 1])
  surv_0 <- sum(weight[data$A == 0] * pseudo[data$A == 0])/sum(weight[data$A == 0])
  spce <- surv_1 - surv_0 
  
  return(spce)
}


####################################################################
##### estimate the SPCE of MiPS-1010 estimator at time point 5 #####
####################################################################
spce.mips1010 <- SPCE.MiPS(data.sim=data,
                           ps1.cov = c('Z1','Z2','Z3'), ps2.cov = NULL,
                           or1.cov = c('Z1','Z2'), or2.cov = NULL,
                           est_time = 5,##the prespecified time point
                           hidden.neurons=c(4,4), learning.rate.global=0.001, momentum.global=0.5)
spce.mips1010##0.1687227
