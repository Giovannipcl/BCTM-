te.BD_pya <- function(y=data[,y],x=data[,x], knotsy, knotsx,  q1, q2, center=T){
  X1 <- splineDesign(knotsy, y, ord = 4, derivs=1, outer.ok=T)
  X2 <- splineDesign(knotsx, x, ord =4, outer.ok=T)
  n <- length(y)
  X <- matrix(0, n, q1 * q2)
  for (i in 1:n) {
    X[i, ] <- X1[i, ] %x% X2[i, ]
  }
  IS <- matrix(0, q1, q1)
  IS[1:q1, 1] <- 1
  for (j in 2:q1) IS[j, 2:j] <- 1
  I <- diag(q2)
  Sig <- IS %x% I
  X <- X %*% Sig
  D <- diag(q1 * q2)
  D <- D[, -q2]
  D1 <- t(diff(diag(q2)))
  D[1:q2, 1:(q2 - 1)] <- D1
  
  if(center) X <- X %*% D
  return(X)
}


MH = function(index, pila, maximos, sigma_ini, lista, lpost, weight){
  dens_norm = marginal(index,pila,maximos,sigma_ini,lista,lpost,weight)
  #dens_norm = marginal(80,post_multi,maximos,sigma_ini,lista,lpost,weight)
  norm = dens_norm$pila_marginal
  xis = dens_norm$xis
  x = rep(0,5000)
  x[1] = xis[norm == max(norm)]     #initialize; I've set arbitrarily set this to 3
  for(i in 2:5000){
    current_x = x[i-1]
    proposed_x = current_x + rnorm(1,mean=0,sd=1)
    #proposed_x = runif(1,current_x -0.2,current_x+0.2) 
    A = min(1,dens_est(proposed_x,norm,xis)/dens_est(current_x,norm,xis))
    #A = interp1(x = xis,y = norm, xi = proposed_x,method = "spline")/interp1(x = xis,y = norm, xi = current_x,method = "spline")
    if(runif(1)<A){
      x[i] = proposed_x       # accept move with probabily min(1,A)
    } else {
      x[i] = current_x        # otherwise "reject" move, and stay where we are
    }
  }
  
  return(x[seq(2000,5000,2)])
}


dens_est = function(x,norm,xis){
  if(x > max(xis) || x< min(xis)){
    return(0)
  }else{
    ind1 = which(abs(xis - x) == sort((abs(xis - x)))[1])
    ind2 = which(abs(xis - x) == sort((abs(xis - x)))[2])
    pesos = abs(x - xis[c(ind1,ind2)])/(sum(abs(x - xis[c(ind1,ind2)])))
    return(interp1(x = xis, y = norm,xi = x,method = "nearest"))}
}


#################WAIC e DIC ################################

ll_dev = function(beta, A,  Al,index){
  gamma = gamma_f(beta,index)
  
  
  res = sum(dnorm(A%*%gamma,0,1,log = TRUE)) + sum(log(Al%*%gamma)) 
  return(-2*res)
  
}

ll_exp = function(beta, A,  Al,index){
  gamma = gamma_f(beta,index)
  
  
  res = (dnorm(A%*%gamma,0,1,log = TRUE)) + (log(Al%*%gamma))
  
  return(exp(res))
  
}

ll2 = function(beta, A,  Al,index){
  gamma = gamma_f(beta,index)
  
  
  res = (dnorm(A%*%gamma,0,1,log = TRUE)) + (log(Al%*%gamma))
  
  return(res)
  
}

ll_2 = function(beta, A,  Al,index){
  gamma = gamma_f(beta,index)
  
  
  res = (dnorm(A%*%gamma,0,1,log = TRUE)) + (log(Al%*%gamma))
  
  return(res^2)
  
}

gamma_f = function(beta,index){
  gammat = vector(length = length(index))
  gammat[index] = exp(beta[index])
  gammat[!index] = beta[!index]
  return(gammat)
}


