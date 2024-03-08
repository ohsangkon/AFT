######################################################################################
################## simulation ###########################################
######################################################################################
library(tidyverse)
library(sn)
library(nnls)
library(latex2exp)
library(e1071)
library(survival)
library(aftgee)
library(ggplot2)
library(emplik)

options(scipen = 100)
######################################################################################
iter = 200
n_sample = 200
beta_true = c(2, 1, -1) 
tau = 1.5

######################c################################################################
## Assumption1. logT ~ Skewt(xi = 0, sigma = 1, lambda = -15, nu = 3)
######################################################################################
#### Generate Data #################################
###################################################
X = list()
logT = list()  # True log survival time
logY = list()  # Observed log survival time
delta = list() # delta = 1: observed, delta = 0: censored
Eps_true = list()

set.seed(1234)
for(i in 1 : iter){

	# Covariate
		X_1 = rnorm(n_sample, 0, 1)
		X_2 = rbinom(n_sample, 1, 0.5)
		X[[i]] = cbind(X_1, X_2)
	
	#True logT
		Eps_true[[i]] = rst(n_sample, xi = 0, omega = 1, alpha = -15, nu = 3)
		Eps_true[[i]] = (Eps_true[[i]] - mean(Eps_true[[i]]))/sd(Eps_true[[i]])

	   design.X = cbind(1, X[[i]])
		logT[[i]] = (design.X %*% beta_true) + Eps_true[[i]]

	#Censor: delta
		logC <- runif(n_sample, 0, tau)	

		delta[[i]] = (logT[[i]] < logC) * 1	

	# censored data 
		logY[[i]] = pmin(logT[[i]], logC)  #logT[[i]]

		if(delta[[i]][which.max(Eps_true[[i]]),] == 0){
			logY[[i]][which.max(Eps_true[[i]]),] = logT[[i]][which.max(Eps_true[[i]]),]
			delta[[i]][which.max(Eps_true[[i]]),] = 1
		}
		
}

## censored rate:
censored_rate = NULL
for(i in 1 : iter){
	censored_rate[i] = mean(delta[[i]] == 0)
}
summary(censored_rate)   

####################################
### Results #######################
###################################

normal = list()     ; beta_normal = list()
skewnormal = list() ; beta_skewnormal = list()
NSNSM = list()      ;  beta_NSNSM = list()

iter = 200
for(i in 1:iter){

	x = X[[i]]
	y = logY[[i]]
	d = delta[[i]]

	normal[[i]] = survreg(Surv(time = exp(y), event = d) ~ x, dist = 'lognormal')
	skewnormal[[i]] = AFT_SN(x = x, y = y, delta = d)
	NSNSM[[i]] = AFT_NSNSM(x, y, d, simulation_index = i)

	beta_normal[[i]]  = normal[[i]]$coefficients %>% unname()
	beta_skewnormal[[i]]  = skewnormal[[i]][[1]]
	beta_NSNSM[[i]]  = NSNSM[[i]][[1]]

	print(i)
}


### Gehan ###
gehan = list()
beta_gehan = list()

for(i in 1 : iter){
	x = X[[i]]
	y = logY[[i]]
	d = delta[[i]]

	dt = as.data.frame(cbind(exp(y), x, d))
	colnames(dt) = c("Y", "X1", "X2", "delta")
	gehan[[i]] =  aftsrr(Surv(Y, delta) ~ X1 + X2, data = dt)
	beta_gehan[[i]] = gehan[[i]]$beta

	esp = y - x %*% as.matrix(gehan[[i]]$beta)
	
	KME <- survfit(Surv(esp, d)~ 1, conf.type="none")
	intercept = KME$time%*%diff(c(0,(1-KME$surv))) # Estimated intercept for Gehan

	beta_gehan[[i]] = c(intercept, beta_gehan[[i]])

	print(i)
}


## GEE ####
gee = list()
beta_gee = list()

for(i in 1 : iter){
	x = X[[i]]
	y = logY[[i]]
	d = delta[[i]]

	dt = as.data.frame(cbind(exp(y), x, d))
	colnames(dt) = c("Y", "X1", "X2", "delta")
	gee[[i]] =  aftgee(Surv(Y, delta) ~ X1 + X2, data = dt)
	beta_gee[[i]] = gee[[i]]$coef.res

	print(i)
}

######################################
#MSE
beta_hat = matrix(0, nrow = iter, ncol = length(beta_true))
L2_normal = beta_hat
L2_skewnormal = beta_hat
L2_NSNSM = beta_hat
L2_gee = beta_hat
L2_gehan = beta_hat
 
for(i in 1 : iter){
	L2_normal[i,] = ((beta_true - beta_normal[[i]])^2)
	L2_skewnormal[i,] = ((beta_true - beta_skewnormal[[i]])^2)
	L2_NSNSM[i,] = ((beta_true - beta_NSNSM[[i]])^2)
	L2_gee[i,] = ((beta_true - beta_gee[[i]])^2)
	L2_gehan[i,] = ((beta_true - beta_gehan[[i]])^2)
}

apply(L2_normal, 2, mean) #* 100
apply(L2_skewnormal, 2, mean) #* 100
apply(L2_NSNSM, 2, mean) #* 100
apply(L2_gee, 2, mean) #* 100
apply(L2_gehan, 2, mean) #* 100

#Bias
bias_normal = beta_hat
bias_skewnormal = beta_hat
bias_NSNSM = beta_hat
bias_gee = beta_hat
bias_gehan = beta_hat

for(i in 1 : iter){
	bias_normal[i,] = beta_normal[[i]]
	bias_skewnormal[i,] = beta_skewnormal[[i]]
	bias_NSNSM[i,] = beta_NSNSM[[i]]
	bias_gee[i,] = beta_gee[[i]]
	bias_gehan[i,] = beta_gehan[[i]]
}

(apply(bias_normal,2,mean) - beta_true) 
(apply(bias_skewnormal,2,mean) - beta_true)
(apply(bias_NSNSM,2,mean) - beta_true) 
(apply(bias_gee,2,mean) - beta_true) 
(apply(bias_gehan,2,mean) - beta_true) 

###########################################################
######################################################################################
## Assumption3. logT ~ Normal(mean = 0, sd = 1)
######################################################################################
#### Generate Data #################################
###################################################
X3 = list()
logT3 = list()  # True log survival time
logY3 = list()  # Observed log survival time
delta3 = list() # delta = 1: observed, delta = 0: censored
Eps_true3 = list()

set.seed(530)
for(i in 1 : iter){

	# Covariate
		X_1 = rnorm(n_sample, 0, 1)
		X_2 = rbinom(n_sample, 1, 0.5)
		X3[[i]] = cbind(X_1, X_2)
		
	#True logT
		Eps_true3[[i]] = rnorm(n_sample, mean = 0, sd = 1)

	   design.X = cbind(1, X3[[i]])
		logT3[[i]] = (design.X %*% beta_true) + Eps_true3[[i]]

	#Censor: delta
		logC <- runif(n_sample, 0, tau)	

		delta3[[i]] = (logT3[[i]] < logC) * 1	

	# censored data 
		logY3[[i]] = pmin(logT3[[i]], logC)  #logT3[[i]]

		if(delta3[[i]][which.max(Eps_true3[[i]]),] == 0){
			logY3[[i]][which.max(Eps_true3[[i]]),] = logT3[[i]][which.max(Eps_true3[[i]]),]
			delta3[[i]][which.max(Eps_true3[[i]]),] = 1
		}
}

## censored rate:
censored_rate = NULL
for(i in 1 : iter){
	censored_rate[i] = mean(delta3[[i]] == 0)
}
summary(censored_rate)   


####################################
### Results #######################
###################################
normal3 = list()     ; beta_normal3 = list()
skewnormal3 = list() ; beta_skewnormal3 = list()
NSNSM3 = list()      ;  beta_NSNSM3 = list()

iter = 200
for(i in 1 : iter){

	x = X3[[i]]
	y = logY3[[i]]
	d = delta3[[i]]

	normal3[[i]] = survreg(Surv(time = exp(y), event = d) ~ x, dist = 'lognormal')
	skewnormal3[[i]] = AFT_SN(x = x, y = y, delta = d)
	NSNSM3[[i]] = AFT_NSNSM(x, y, d, simulation_index = i)

	beta_normal3[[i]]  = normal3[[i]]$coefficients %>% unname()
	beta_skewnormal3[[i]]  = skewnormal3[[i]][[1]]
	beta_NSNSM3[[i]]  = NSNSM3[[i]][[1]]

	print(i)
}

### Gehan ###
gehan3 = list()
beta_gehan3 = list()

for(i in 1 : iter){
	x = X3[[i]]
	y = logY3[[i]]
	d = delta3[[i]]

	dt = as.data.frame(cbind(exp(y), x, d))
	colnames(dt) = c("Y", "X1", "X2", "delta")
	gehan3[[i]] =  aftsrr(Surv(Y, delta) ~ X1 + X2, data = dt)
	beta_gehan3[[i]] = gehan3[[i]]$beta

	esp = y - x %*% as.matrix(gehan3[[i]]$beta)
	
	KME <- survfit(Surv(esp, d)~ 1, conf.type="none")
	intercept = KME$time%*%diff(c(0,(1-KME$surv))) # Estimated intercept for Gehan

	beta_gehan3[[i]] = c(intercept, beta_gehan3[[i]])

	print(i)
}


## GEE ####
gee3 = list()
beta_gee3 = list()

for(i in 1 : iter){
	x = X3[[i]]
	y = logY3[[i]]
	d = delta3[[i]]

	dt = as.data.frame(cbind(exp(y), x, d))
	colnames(dt) = c("Y", "X1", "X2", "delta")
	gee3[[i]] =  aftgee(Surv(Y, delta) ~ X1 + X2, data = dt)
	beta_gee3[[i]] = gee3[[i]]$coef.res

	print(i)
}

###################################################
#MSE
beta_hat = matrix(0, nrow = iter, ncol = length(beta_true))
L2_normal3 = beta_hat
L2_skewnormal3 = beta_hat
L2_NSNSM3 = beta_hat
L2_gee3 = beta_hat
L2_gehan3 = beta_hat
 
for(i in 1 : iter){
	L2_normal3[i,] = ((beta_true - beta_normal3[[i]])^2)
	L2_skewnormal3[i,] = ((beta_true - beta_skewnormal3[[i]])^2)
	L2_NSNSM3[i,] = ((beta_true - beta_NSNSM3[[i]])^2)
	L2_gee3[i,] = ((beta_true - beta_gee3[[i]])^2)
	L2_gehan3[i,] = ((beta_true - beta_gehan3[[i]])^2)
}

apply(L2_normal3, 2, mean) #* 100
apply(L2_skewnormal3, 2, mean)#* 100
apply(L2_NSNSM3, 2, mean)#* 100
apply(L2_gee3, 2, mean)#* 100
apply(L2_gehan3, 2, mean)#* 100

#Bias
bias_normal3 = beta_hat
bias_skewnormal3 = beta_hat
bias_NSNSM3 = beta_hat
bias_gee3 = beta_hat
bias_gehan3 = beta_hat

for(i in 1 : iter){
	bias_normal3[i,] = beta_normal3[[i]]
	bias_skewnormal3[i,] = beta_skewnormal3[[i]]
	bias_NSNSM3[i,] = beta_NSNSM3[[i]]
	bias_gee3[i,] = beta_gee3[[i]]
	bias_gehan3[i,] = beta_gehan3[[i]]
}

(apply(bias_normal3,2,mean) - beta_true) 
(apply(bias_skewnormal3,2,mean) - beta_true)
(apply(bias_NSNSM3,2,mean) - beta_true)
(apply(bias_gee3,2,mean) - beta_true)
(apply(bias_gehan3,2,mean) - beta_true)

###########################################################
######################################################################################
## Assumption4. logT ~ t(df = 3)
######################################################################################
#### Generate Data #################################
###################################################
X4 = list()
logT4 = list()  # True log survival time
logY4 = list()  # Observed log survival time
delta4 = list() # delta = 1: observed, delta = 0: censored
Eps_true4 = list()

set.seed(86542331)
for(i in 1 : iter){

	# Covariate
		X_1 = rnorm(n_sample, 0, 1)
		X_2 = rbinom(n_sample, 1, 0.5)
		X4[[i]] = cbind(X_1, X_2)
		
	#True logT
		Eps_true4[[i]] = rt(n_sample, df = 3)
		Eps_true4[[i]] = (Eps_true4[[i]])/sd(Eps_true4[[i]])

	   design.X = cbind(1, X4[[i]])
		logT4[[i]] = (design.X %*% beta_true) + Eps_true4[[i]]

	#Censor: delta
		logC <- runif(n_sample, 0, tau)	

		delta4[[i]] = (logT4[[i]] < logC) * 1	

	# censored data 
		logY4[[i]] = pmin(logT4[[i]], logC)  #logT4[[i]]

		if(delta4[[i]][which.max(Eps_true4[[i]]),] == 0){
			logY4[[i]][which.max(Eps_true4[[i]]),] = logT4[[i]][which.max(Eps_true4[[i]]),]
			delta4[[i]][which.max(Eps_true4[[i]]),] = 1
		}

}

## censored rate:
censored_rate = NULL
for(i in 1 : iter){
	censored_rate[i] = mean(delta4[[i]] == 0)
}
summary(censored_rate)   


####################################
### Results #######################
###################################
normal4 = list()     ; beta_normal4 = list()
skewnormal4 = list() ; beta_skewnormal4 = list()
NSNSM4 = list()      ;  beta_NSNSM4 = list()


iter = 200
for(i in 1 : iter){

	x = X4[[i]]
	y = logY4[[i]]
	d = delta4[[i]]

	normal4[[i]] = survreg(Surv(time = exp(y), event = d) ~ x, dist = 'lognormal')
	skewnormal4[[i]] = AFT_SN(x = x, y = y, delta = d)
	NSNSM4[[i]] = AFT_NSNSM(x, y, d, simulation_index = i)

	beta_normal4[[i]]  = normal4[[i]]$coefficients %>% unname()
	beta_skewnormal4[[i]]  = skewnormal4[[i]][[1]]
	beta_NSNSM4[[i]]  = NSNSM4[[i]][[1]]

	print(i)
}

### Gehan ###
gehan4 = list()
beta_gehan4 = list()

for(i in 1 : iter){
	x = X4[[i]]
	y = logY4[[i]]
	d = delta4[[i]]

	dt = as.data.frame(cbind(exp(y), x, d))
	colnames(dt) = c("Y", "X1", "X2", "delta")
	gehan4[[i]] =  aftsrr(Surv(Y, delta) ~ X1 + X2, data = dt)
	beta_gehan4[[i]] = gehan4[[i]]$beta

	esp = y - x %*% as.matrix(gehan4[[i]]$beta)
	
	KME <- survfit(Surv(esp, d)~ 1, conf.type="none")
	intercept = KME$time%*%diff(c(0,(1-KME$surv))) # Estimated intercept for Gehan

	beta_gehan4[[i]] = c(intercept, beta_gehan4[[i]])

	print(i)
}



## GEE ####
gee4 = list()
beta_gee4 = list()

for(i in 1 : iter){
	x = X4[[i]]
	y = logY4[[i]]
	d = delta4[[i]]

	dt = as.data.frame(cbind(exp(y), x, d))
	colnames(dt) = c("Y", "X1", "X2", "delta")
	gee4[[i]] =  aftgee(Surv(Y, delta) ~ X1 + X2, data = dt)
	beta_gee4[[i]] = gee4[[i]]$coef.res

	print(i)
}

#######################################################
#MSE
beta_hat = matrix(0, nrow = iter, ncol = length(beta_true))
L2_normal4 = beta_hat
L2_skewnormal4 = beta_hat
L2_NSNSM4 = beta_hat
L2_gee4 = beta_hat
L2_gehan4 = beta_hat

for(i in 1 : iter){
	L2_normal4[i,] = ((beta_true - beta_normal4[[i]])^2)
	L2_skewnormal4[i,] = ((beta_true - beta_skewnormal4[[i]])^2)
	L2_NSNSM4[i,] = ((beta_true - beta_NSNSM4[[i]])^2)
	L2_gee4[i,] = ((beta_true - beta_gee4[[i]])^2)
	L2_gehan4[i,] = ((beta_true - beta_gehan4[[i]])^2)
}

apply(L2_normal4, 2, mean)
apply(L2_skewnormal4, 2, mean)
apply(L2_NSNSM4, 2, mean) 
apply(L2_gee4, 2, mean) 
apply(L2_gehan4, 2, mean) 

#Bias
bias_normal4 = beta_hat
bias_skewnormal4 = beta_hat
bias_NSNSM4 = beta_hat
bias_gee4 = beta_hat
bias_gehan4 = beta_hat

for(i in 1 : iter){
	bias_normal4[i,] = beta_normal4[[i]]
	bias_skewnormal4[i,] = beta_skewnormal4[[i]]
	bias_NSNSM4[i,] = beta_NSNSM4[[i]]
	bias_gee4[i,] = beta_gee4[[i]]
	bias_gehan4[i,] = beta_gehan4[[i]]
}

(apply(bias_normal4,2,mean) - beta_true) 
(apply(bias_skewnormal4,2,mean) - beta_true)
(apply(bias_NSNSM4,2,mean) - beta_true)
(apply(bias_gee4,2,mean) - beta_true) 
(apply(bias_gehan4,2,mean) - beta_true) 

######################################################################################
## Assumption5. logT ~ Gumbel(0,5)
######################################################################################
#### Generate Data #################################
###################################################
X5 = list()
logT5 = list()  # True log survival time
logY5 = list()  # Observed log survival time
delta5 = list() # delta = 1: observed, delta = 0: censored
Eps_true5 = list()

set.seed(5302)
for(i in 1 : iter){
	library(evd)
	# Covariate
		X_1 = rnorm(n_sample, 0, 1)
		X_2 = rbinom(n_sample, 1, 0.5)
		X5[[i]] = cbind(X_1, X_2)
		
	#True logT
		Eps_true5[[i]] = rgumbel(n_sample, loc = 0, scale = 5)
		Eps_true5[[i]] = (Eps_true5[[i]] - mean(Eps_true5[[i]]))/sd(Eps_true5[[i]])

	   design.X = cbind(1, X5[[i]])
		logT5[[i]] = (design.X %*% beta_true) + Eps_true5[[i]]

	#Censor: delta
		logC <- runif(n_sample, 0, tau)	

		delta5[[i]] = (logT5[[i]] < logC) * 1	

	# censored data 
		logY5[[i]] = pmin(logT5[[i]], logC)  #logT5[[i]]

		if(delta5[[i]][which.max(Eps_true5[[i]]),] == 0){
			logY5[[i]][which.max(Eps_true5[[i]]),] = logT5[[i]][which.max(Eps_true5[[i]]),]
			delta5[[i]][which.max(Eps_true5[[i]]),] = 1
		}

}

## censored rate:
censored_rate = NULL
for(i in 1 : iter){
	censored_rate[i] = mean(delta5[[i]] == 0)
}
summary(censored_rate)   


####################################
### Results #######################
###################################
normal5 = list()     ; beta_normal5 = list()
skewnormal5 = list() ; beta_skewnormal5 = list()
NSNSM5 = list()      ;  beta_NSNSM5 = list()

iter = 200
for(i in 1: iter){

	x = X5[[i]]
	y = logY5[[i]]
	d = delta5[[i]]

	normal5[[i]] = survreg(Surv(time = exp(y), event = d) ~ x, dist = 'lognormal')
	skewnormal5[[i]] = AFT_SN(x = x, y = y, delta = d)
	NSNSM5[[i]] = AFT_NSNSM(x, y, d, simulation_index = i)

	beta_normal5[[i]]  = normal5[[i]]$coefficients %>% unname()
	beta_skewnormal5[[i]]  = skewnormal5[[i]][[1]]
	beta_NSNSM5[[i]]  = NSNSM5[[i]][[1]]

	print(i)
}


### Gehan ###
gehan5 = list()
beta_gehan5 = list()

for(i in 1 : iter){
	x = X5[[i]]
	y = logY5[[i]]
	d = delta5[[i]]

	dt = as.data.frame(cbind(exp(y), x, d))
	colnames(dt) = c("Y", "X1", "X2", "delta")
	gehan5[[i]] =  aftsrr(Surv(Y, delta) ~ X1 + X2, data = dt)
	beta_gehan5[[i]] = gehan5[[i]]$beta

	esp = y - x %*% as.matrix(gehan5[[i]]$beta)
	
	KME <- survfit(Surv(esp, d)~ 1, conf.type="none")
	intercept = KME$time%*%diff(c(0,(1-KME$surv))) # Estimated intercept for Gehan

	beta_gehan5[[i]] = c(intercept, beta_gehan5[[i]])

	print(i)
}



## GEE ####
gee5 = list()
beta_gee5 = list()

for(i in 1 : iter){
	x = X5[[i]]
	y = logY5[[i]]
	d = delta5[[i]]

	dt = as.data.frame(cbind(exp(y), x, d))
	colnames(dt) = c("Y", "X1", "X2", "delta")
	gee5[[i]] =  aftgee(Surv(Y, delta) ~ X1 + X2, data = dt)
	beta_gee5[[i]] = gee5[[i]]$coef.res

	print(i)
}

#######################################################
#MSE
beta_hat = matrix(0, nrow = iter, ncol = length(beta_true))
L2_normal5 = beta_hat
L2_skewnormal5 = beta_hat
L2_NSNSM5 = beta_hat
L2_gee5 = beta_hat
L2_gehan5 = beta_hat

for(i in 1 : iter){
	L2_normal5[i,] = ((beta_true - beta_normal5[[i]])^2)
	L2_skewnormal5[i,] = ((beta_true - beta_skewnormal5[[i]])^2)
	L2_NSNSM5[i,] = ((beta_true - beta_NSNSM5[[i]])^2)
	L2_gee5[i,] = ((beta_true - beta_gee5[[i]])^2)
	L2_gehan5[i,] = ((beta_true - beta_gehan5[[i]])^2)
}

apply(L2_normal5, 2, mean) 
apply(L2_skewnormal5, 2, mean) 
apply(L2_NSNSM5, 2, mean) 
apply(L2_gee5, 2, mean) 
apply(L2_gehan5, 2, mean) 

#Bias
bias_normal5 = beta_hat
bias_skewnormal5 = beta_hat
bias_NSNSM5 = beta_hat
bias_gee5 = beta_hat
bias_gehan5 = beta_hat

for(i in 1 : iter){
	bias_normal5[i,] = beta_normal5[[i]]
	bias_skewnormal5[i,] = beta_skewnormal5[[i]]
	bias_NSNSM5[i,] = beta_NSNSM5[[i]]
	bias_gee5[i,] = beta_gee5[[i]]
	bias_gehan5[i,] = beta_gehan5[[i]]

}

(apply(bias_normal5,2,mean) - beta_true) 
(apply(bias_skewnormal5,2,mean) - beta_true) 
(apply(bias_NSNSM5,2,mean) - beta_true) 
(apply(bias_gee5,2,mean) - beta_true) 
(apply(bias_gehan5,2,mean) - beta_true) 

###########################################################
