################################3
### Loading the packages
library(tidyverse)
library(dplyr)
library(Matrix)
library(splines)
library(splines2)
library(MASS)
library(mgcv)
library(scam)
library(scales)
library(invgamma)
library(pracma)

dir <- dirname(rstudioapi::getActiveDocumentContext()$path)  
dir_models <- file.path(dir, "models")  
dir_cpp <- file.path(dir, "cpp")
### Loading the original matrix model to ensure that we are using the same data
# The fram_te_m10_re is originally from "Carlan et al. [2024], Bayesian Conditional Transformation Models"
setwd(dir_models)
load("fram_te_m10.RData") # processed data from Carlan. It is avaiable at " https://github.com/manucarl/BCTM/tree/main/processed_data/framingham"



## df and db are the original data 
df <- read.csv("https://raw.githubusercontent.com/manucarl/BCTM/main/data/framingham.csv")
head(df)
attach(df)
db <- read.csv("https://raw.githubusercontent.com/manucarl/BCTM/main/data/framingham.csv")
##################################


########### Getting auxiliaries functions and loading the .cpp codes
source("All_models.R")
source("funcs_aux.R")

library(Rcpp)
library(RcppArmadillo)
library(RcppParallel)
setwd(dir_cpp)
sourceCpp("model2.cpp")


#################################################################################
################################
## Rescale the data
df$age <- rescale(df$age, to = c(0, 1))
df$year <- rescale(df$year, to = c(0, 1))
#standardize
df$cholst <- scale(df$cholst)



q = c(10,10)  # size of y (cholst) and x (age) basis
X = model.matrix(~1+sex + year,data = df) #fixed effects
y = df$cholst # response variable
head(df)

### Crearing the tensor between cholesterol and age (this code is optinal since we are getting all these matrices from
                                                      #Carlan et al. [2024], Bayesian Conditional Transformation Models)
sm <- smoothCon(s(cholst,age,bs="tesmi1",k=q),data=df,
                scale.penalty=F, absorb.cons=TRUE, knots=NULL)[[1]]

# Getting derivative of the transformation funcion
Bcl = te.BD_pya(y,df$age,sm$knots[[1]],sm$knots[[2]], q[1],q[2])
Xi0 = matrix(0, ncol = ncol(X),nrow = nrow(X))


## matrix of the transformation functions from Carlan et al. [2024], Bayesian Conditional Transformation Models, and fixed effects X
A = cbind(object_te$X[,2:100],X)
Ad = cbind(object_te$Xp[,2:100],Xi0)  # derivative of the transformation functions


## Penalization matrices for the prior distributions
K1 = sm$S[[1]]
K2 = sm$S[[2]]


## Index of the parameters (TRUE for the exponentied parameters, and FALSE otherwise)
index1 = sm$p.ident
indext = c(index1,rep(FALSE,ncol(X)))

# auxiliar matrix
K3 = matrix(0,2,2)


## Fitting the model ##################################################
##################################################
system.time({m12 = bctm_tensor_multi(A,Ad,X,K1,K2,K3,q,indext)})
##################################################
##################################################


#### Auxiliar functions ##################################################
marginal = function(index,pila,maximos,sigma_ini,lista,lpost,weight){
  pila_marginal = vector()
  xis = xis_c2(index,maximos,sigma_ini)
  mar = function(j){
    return(sum(sapply(lapply(1:(length(lista)), function(i) pila[[i]][,index]),'[[',j)*(lpost*weight)))
  }
  pila_marginal = unlist(map(1:length(xis), mar))
  return(list(pila_marginal = as.vector(pila_marginal), xis = xis))
}

xis_c2 = function(index,max,sigma_ini){
  if(index %in% which(indext == FALSE)){
    return(seq(max[index]-2.5*sqrt(sigma_ini[index]),max[index]+2.5*sqrt(sigma_ini[index]),length = 31))
  }else{
    return(seq(max[index]-2.5*sqrt(sigma_ini[index]),max[index]+2.5*sqrt(sigma_ini[index]),length = 31))
  }
}


#################################################################

post_multi = list()
points_t = m12[[1]]$sigma   #extracting the points explored from the hyperparameter posterior distribution
maximos = m12[[2]]$maximos  # getting the maximum of the parameters posterior pdf 
sigma_ini = m12[[2]]$sigma_ini  # Diagonal of the Hessian amtrix of the laplace approximation for (beta|y, tau), 
#it is the variance parameters of the approximated distribution
lista = as.list(as.data.frame(t(points_t))) 


unl = (sapply(m12, '[[', 'unn_log_post')) - mean(sapply(m12, `[[`, 'unn_log_post'))

weight = c(abs(diff(points_t[,2])[1])*abs(diff(points_t[,1],7)[1])) # weight for the integral of the hyperparameter distributions
post_alpha  = exp(unl)/(sum(exp(unl)*weight)) ### Hyperparameter posterior distribution normalized
lpost = post_alpha 

lpost2 <- cut((post_alpha), 10, label = FALSE)
colr <- rev(heat.colors(10))

### Plot of the hyperparameter distribution
plot(points_t[,1] ~ points_t[,2],
     col = colr[lpost2], pch = 20)

pontos = points_t

lpost2_norm = lpost

### extracting the Laplace approximation for (beta|y, tau)
for(j in 1:length(pontos[,1])){
  post_multi[[j]] = matrix(unlist(m12[[j]]$posti_cond),nrow = 31,ncol = ncol(A))
}


### Getting the marginals of the posterior parameters
marginals2 = lapply(1:length(maximos),function(ii) marginal(ii,post_multi, maximos = maximos,
                                                           sigma_ini = sigma_ini,lista = lista,lpost = lpost,weight = weight))



### Simulating values from these marginals
Beta = mclapply(seq(1,length(maximos),1),MH,post_multi,maximos,sigma_ini,lista,lpost,weight,mc.cores =1 )
Beta = matrix(unlist(Beta),ncol = ncol(A),nrow = length(Beta[[1]]))



################ WAIC and DIC ##############################
i = 23 #(49/2), mode of hyperparameter
mu_gaus = m12[[i]]$x0
sigma_gaus = solve(m12[[i]]$H1)
sam = mvrnorm(1000,mu_gaus,(sigma_gaus))

py = vector()
Tj = matrix(0,length(y),length(m12))
Tj2 = matrix(0,length(y),length(m12))
inte = vector()
inte2 = vector()
set.seed(1)

for(i in 1:length(m12)){
  print(i)
  mu_gaus = m12[[i]]$x0
  sigma_gaus = (m12[[i]]$H1)
  chol_H_ = chol(sigma_gaus)
  sigma_gaus = chol2inv(chol_H_)
  sam = mvrnorm(1000,mu_gaus,(sigma_gaus))
  Tj[,i] = apply((apply((sam),1,ll_exp, A = A, Al = Ad,index = indext)),1,mean)
  Tj2[,i] = apply(apply(sam,1,ll_2,A,Ad,indext),1,mean) - apply(apply(sam,1,ll2,A,Ad,indext),1,mean)^2
  inte[i] = mean(apply((sam),1,ll_dev, A = A, Al = Ad,index = indext))
  inte2[i] = ll_dev(mu_gaus,A,Ad,indext)
}

WAIC = -2*sum(log(Tj%*%(lpost2_norm*weight))) + 2*sum((Tj2%*%(lpost2_norm*weight)));WAIC

inte_int = sum(inte*lpost2_norm*weight)
beta_mean = apply(Beta,2,mean)
inte2 = ll_dev(beta_mean,A,Ad, indext)
pd = inte_int - inte2;pd
DIC =inte2  + pd;DIC

#####################################################

Gamma = apply(Beta,1,gamma_f,index = indext)

# function to estimate the conditional distribution functions
prob_estim = function(ind,Gamma){
  gamma = t(Gamma)
  Z = sm$Zc
  ynew_o = seq(min(db$cholst)-10,max(db$cholst)+10,1)
  
  ynew = (ynew_o - attr(y,"scaled:center"))/(attr(y,"scaled:scale"))

  
  dfnew = as.data.frame(cbind(ynew, df[ind,6]))
  names(dfnew) = c("cholst","age")
  B = PredictMat(sm, as.data.frame(dfnew))%*%Z
  Bl = (te.BD_pya(dfnew$cholst,dfnew$age,sm$knots[[1]],sm$knots[[2]], q[1],q[2]))
  dfnew = as.data.frame(cbind(ynew, df[ind,6]))
  names(dfnew) = c("cholst","age")
  
  bt_samples <-  object_te$samples$beta[2000:4000,]
  bt_samples[,object_te$model$exp_ident] <- exp(bt_samples[,object_te$model$exp_ident])
  
  h1 = B%*%(t(gamma[,1:(q[1]*q[2] - 1)])) 
  h1c = B%*%(t(bt_samples[,2:(q[1]*q[2])])) 
  h1l = Bl%*%(t(gamma[,1:(q[1]*q[2] - 1)])) 
  h1lc = Bl%*%(t(bt_samples[,2:(q[1]*q[2])])) 
  
  
  h2 = X[ind,]%*%t(gamma[,(q[1]*q[2]):ncol(gamma)])
  h2c = X[ind,]%*%t(bt_samples[,c(1,(q[1]*q[2] + 1):ncol(gamma))])
  h1 = apply(h1,1,mean)
  h1c = apply(h1c,1,mean)
  h1l = apply(h1l,1,mean)
  h1lc = apply(h1lc,1,mean)
  h2 = mean(h2)
  h2c = mean(h2c)
  
  Fa = pnorm(h1 + h2)
  Fac = pnorm(h1c + h2c)
  #plot(Fa,type = "l")
  #points(Fac, type = "l",col = "blue")
  Faord = sort(Fa)
  Facord = sort(Fac)
  qi = (ynew_o[Faord > 0.05])[1]
  qs = (ynew_o[Faord > 0.95])[1]
  qm = (ynew_o[Faord > 0.5])[1]
  qic = (ynew_o[Facord > 0.05])[1]
  qsc = (ynew_o[Facord > 0.95])[1]
  qmc = (ynew_o[Facord > 0.5])[1]

  
  return(list(dif = abs(Fa- Fac), qi = qi, qic = qic, qs = qs, qsc =qsc,qm = qm, qmc = qmc))
  
  
}


Estim = lapply(seq(1:(dim(df)[1])),prob_estim, Gamma = Gamma)


qm = sapply(Estim, "[[","qm")     # Median
qmc = sapply(Estim, "[[","qmc")   # Median of original paper
qi = sapply(Estim, "[[","qi")     # 0.05 quantile
qs = sapply(Estim, "[[","qs")     # 0.95 quantile
qic = sapply(Estim, "[[","qic")   # 0.05 quantile of original paper  
qsc = sapply(Estim, "[[","qsc")   #0.95 quantile of original paper
dif = sapply(Estim, '[[',"dif")   # dif between estimated functions distributions from both approaches

dif_df = as.data.frame(cbind(apply(dif, 2, mean),apply(dif, 2, min),apply(dif, 2, max)))
names(dif_df) = c("Mean", "Minimum","Maximum")

library("reshape2") 
library("ggplot2")
data_long <- melt(dif_df)       # Reshaping data frame
head(data_long)


## This plot show the distributions of the minimum, average, and maximum differences from the estimated distributions
ggplot(data_long, aes(variable,value))+
  geom_boxplot(fill = "lightblue", colour = "black")+
  theme_minimal()+ scale_x_discrete(name="") +
  scale_y_continuous(name="Difference", lim = c(0,0.2))+
  theme(axis.title.x = element_text(size = 13),
        axis.title.y= element_text(size = 13),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13))


datapred = as.data.frame(cbind(db$cholst,qmc,qm,qi,qs,qic,qsc))
cr = 1 - sum(datapred$V1 > datapred$qs | datapred$V1 < datapred$qi)/(nrow(datapred))

crc = 1 - sum(datapred$V1 > datapred$qsc | datapred$V1 < datapred$qic)/(nrow(datapred))
# Coverage rate from our approach (cr) and from original paper (crc)
cr;crc


head(datapred)
datapred = datapred[order(datapred$qmc),]
plot(datapred$qmc,type = "l", col =  "darkred",lwd = 2,ylim = c(100,450),lty =2, cex.axis = 1.5,cex.lab = 1.5, ylab = "Cholesterol level")
points(datapred$qsc, type = "l",col = "darkred",lwd = 2)
points(datapred$qic, type = "l",col = "darkred",lwd = 2)
points(datapred$V1,pch = 20)


datapred = as.data.frame(cbind(db$cholst,qmc,qm,qi,qs,qic,qsc))
datapred = datapred[order(datapred$qm),]
plot(datapred$qm,type = "l", col =  "darkred",ylim = c(100,450),lwd = 2,lty =2, cex.axis = 1.5,cex.lab = 1.5, ylab = "Cholesterol level")
points(datapred$qm, type = "l",col = "darkblue",lwd = 2)
points(datapred$qs, type = "l",col = "darkblue",lwd = 2)
points(datapred$qi, type = "l",col = "darkblue",lwd = 2)
points(datapred$V1,pch = 20)

