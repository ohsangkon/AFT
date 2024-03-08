###########################################################
## (Function) AFT with skew normal error  ##
###########################################################
###  Main function
#############################################
library(tidyverse)
library(sn)
library(e1071)
library(survival)

AFT_SN = function(x, y, delta = NULL){

	##initial value setting: 
	if(is.null(delta) == TRUE){
		delta = matrix(1, nrow = nrow(x), ncol = 1)
	}

	# beta : OLS
	lognormal = survreg(Surv(time = exp(y), event = delta) ~ x, dist = 'lognormal')
	beta_hat  = lognormal$coefficients %>% unname()

	# sigma & lambda: method of moment estimator (MME)
	x_obs = x[which(delta == 1) ,]  #observed covariate
	y_obs = y[which(delta == 1) ,]  #observed survival time
	X_obs = cbind(1,x_obs)  # design matrix
	n = length(y_obs)
	Eps = y_obs - X_obs %*% beta_hat
	m2 = sum( (Eps - mean(Eps))^2 ) / (n - 1)
	m3 = sum( (Eps - mean(Eps))^3 ) / (n - 1)
	a = sqrt(2/pi)
	b = (4/pi - 1) * a

	#sigma
	sigma_hat = sqrt(m2 + a^2 * ((m3/b)^(2))^(1/3))

	# lambda
	delta_lambda = ( a^2 + m2 * ((b/m3)^(2))^(1/3) )^(-0.5)
	
	if(skewness(Eps) > 0){
		lambda_hat = sqrt(delta_lambda^2 / (delta_lambda^2 + 1))
	}else{
		lambda_hat = -sqrt(delta_lambda^2 / (delta_lambda^2 + 1))
	}
	Estimates = c(beta_hat, sigma_hat ,lambda_hat)
	Estimates
	
	#Results
	log_likelihood_old_reg = loglikelihood_AFT_SN_Optim(x = x, y = y, delta, Estimates)

	optim_result = 
  		  optim(par = Estimates,
     		     fn = loglikelihood_AFT_SN_Optim,
        		  x = x, y = y, delta = delta, 
     		     method = "BFGS", 
      	     control = list("fnscale" = -1))

	p = ncol(as.matrix(x))   	 
	beta_new = optim_result$par[1:(p+1)]
	sigma_new = optim_result$par[p+2]
	lambda_new = optim_result$par[p+3]
	
	X = cbind(1,x)
	Eps_new <<- y - X %*% beta_new

	Q_support = list()
	Q_support[[1]] = list(prob = 1, xi = 0, sigma = sigma_new, lambda = lambda_new)
	xi_new = -mean_Q_support(Q_support, 0)

	result = list()
	result[[1]] = beta_new
	result[[2]] = sigma_new
	result[[3]] = lambda_new
	result[[4]] = xi_new	
	
	names(result) = c("beta_hat", "sigma_hat", "lambda_hat", "xi_hat")
	return(result)
}


###########################################################
## (Function) Density & Loglikelihood  ##
###########################################################

## 1. Density of AFT_SN
density_AFT_SN = function(Eps, delta, sigma, lambda)
{
	observed = Eps[which(delta == 1)]
	censored = Eps[which(delta != 1)]
	likelihood = rep(0, length(Eps))

  ## E(Eps) = 0 
  ###--->  xi = - sqrt(2/pi) * (lambda / sqrt(1 + lambda^2)) * sigma
  
   xi = - sqrt(2/pi) * (lambda / sqrt(1 + lambda^2)) * sigma

	if(length(observed) != 0 ){
 	  likelihood_observed = dsn(x = observed, 
  				           		    xi = xi,
     					     		    omega = sigma,
        		   		    		 alpha = lambda)

		likelihood_observed[likelihood_observed == 0] = 0.000001
		
		likelihood[which(delta == 1)] = likelihood_observed	
	}

	if(length(censored) != 0){
		likelihood_censored = (1 - psn(x = censored,
                                     xi = xi,
                                     omega = sigma,
                                     alpha = lambda))

		likelihood_censored[likelihood_censored == 0] = 0.000001

		likelihood[which(delta != 1)] = likelihood_censored
   }

  return(likelihood)
}


## 2. loglikelihood 
loglikelihood_AFT_SN = function(Eps, delta, sigma, lambda)
{

  likelihood = density_AFT_SN(Eps = Eps, delta, sigma = sigma, lambda = lambda)
  
  log_likelihood = likelihood %>% log()%>% sum()
  
  return(log_likelihood)
}


## 3. Log_likelihood with Regression
## This is for an objective function for the function 'optim'
# Estimates = c(beta_hat, sigma_hat, lambda_hat)

loglikelihood_AFT_SN_Optim = function(x, y, delta, Estimates)
{
  #x, y: matrix
  p = ncol(as.matrix(x))   	 
  beta_temp = Estimates[1:(p+1)]
  sigma_temp = Estimates[p+2]
  lambda_temp = Estimates[p+3]

  X = cbind(1,x)
  Eps_temp = as.vector(y - X %*% beta_temp)

  log_lik = 
    loglikelihood_AFT_SN(Eps = Eps_temp, delta = delta, sigma = sigma_temp, lambda = lambda_temp)
                        
  return(log_lik)
}



######################################