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
library(moments)
library(emplik)

options(scipen = 100)
######################################################################################
iter = 200
n_sample = 200
beta_true = c(2, 1, -1) 
tau = 1.5

######################c################################################################
## Assumption1. logT ~ Skewt(xi = 0, sigma = 1, lambda = -1, nu = 3)
######################################################################################
#### Generate Data #################################
###################################################
X = list()
logT = list()  # True log survival time
logY = list()  # Observed log survival time
delta = list() # delta = 1: observed, delta = 0: censored
Eps_true = list()

set.seed(123854132)
for(i in 1 : iter){

	# Covariate
		X_1 = rnorm(n_sample, 0, 1)
		X_2 = rbinom(n_sample, 1, 0.5)
		X[[i]] = cbind(X_1, X_2)
	
	#True logT
		Eps_true[[i]] = rst(n_sample, xi = 0, omega = 1, alpha = -1, nu = 3)
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
#L2 norm

L2_normal = NULL 
L2_skewnormal = NULL
L2_gehan = NULL
L2_gee = NULL
L2_NSNSM = NULL


for(i in 1 : iter){
	L2_normal = append(L2_normal, sqrt(sum((beta_normal[[i]] - beta_true)^2)))
	L2_skewnormal = append(L2_skewnormal, sqrt(sum((beta_skewnormal[[i]] - beta_true)^2)))
	L2_gehan = append(L2_gehan, sqrt(sum((beta_gehan[[i]] - beta_true)^2)))
	L2_gee = append(L2_gee, sqrt(sum((beta_gee[[i]] - beta_true)^2)))
	L2_NSNSM = append(L2_NSNSM, sqrt(sum((beta_NSNSM[[i]] - beta_true)^2)))
}


boxplot(L2_normal, L2_skewnormal, L2_gehan, L2_gee, L2_NSNSM, 
		yaxt = "n", col = c(rep("grey",4), rep("pink", 1)), ylim = c(0,2.5))
axis(1, at = 1:5, labels = c("Normal","SN", "Gehan", "GEE", "NSNSM"), cex.axis = 5)
axis(2, cex.axis = 3)


######################################
######################c################################################################
## Assumption2. logT ~ Skewt(xi = 0, sigma = 1, lambda = -10, nu = 3)
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
		Eps_true[[i]] = rst(n_sample, xi = 0, omega = 1, alpha = -10, nu = 3)
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
#L2 norm

L2_normal = NULL 
L2_skewnormal = NULL
L2_gehan = NULL
L2_gee = NULL
L2_NSNSM = NULL


for(i in 1 : iter){
	L2_normal = append(L2_normal, sqrt(sum((beta_normal[[i]] - beta_true)^2)))
	L2_skewnormal = append(L2_skewnormal, sqrt(sum((beta_skewnormal[[i]] - beta_true)^2)))
	L2_gehan = append(L2_gehan, sqrt(sum((beta_gehan[[i]] - beta_true)^2)))
	L2_gee = append(L2_gee, sqrt(sum((beta_gee[[i]] - beta_true)^2)))
	L2_NSNSM = append(L2_NSNSM, sqrt(sum((beta_NSNSM[[i]] - beta_true)^2)))
}


boxplot(L2_normal, L2_skewnormal, L2_gehan, L2_gee, L2_NSNSM, 
		yaxt = "n", col = c(rep("grey",4), rep("pink", 1)), ylim = c(0,2.5))
axis(1, at = 1:5, labels = c("Normal","SN", "Gehan", "GEE", "NSNSM"), cex.axis = 5)
axis(2, cex.axis = 3)


######################################
######################c################################################################
## Assumption3. logT ~ Skewt(xi = 0, sigma = 1, lambda = -50, nu = 3)
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
		Eps_true[[i]] = rst(n_sample, xi = 0, omega = 1, alpha = -50, nu = 3)
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
#L2 norm

L2_normal = NULL 
L2_skewnormal = NULL
L2_gehan = NULL
L2_gee = NULL
L2_NSNSM = NULL


for(i in 1 : iter){
	L2_normal = append(L2_normal, sqrt(sum((beta_normal[[i]] - beta_true)^2)))
	L2_skewnormal = append(L2_skewnormal, sqrt(sum((beta_skewnormal[[i]] - beta_true)^2)))
	L2_gehan = append(L2_gehan, sqrt(sum((beta_gehan[[i]] - beta_true)^2)))
	L2_gee = append(L2_gee, sqrt(sum((beta_gee[[i]] - beta_true)^2)))
	L2_NSNSM = append(L2_NSNSM, sqrt(sum((beta_NSNSM[[i]] - beta_true)^2)))
}


boxplot(L2_normal, L2_skewnormal, L2_gehan, L2_gee, L2_NSNSM, 
		yaxt = "n", col = c(rep("grey",4), rep("pink", 1)), ylim = c(0,2.5))
axis(1, at = 1:5, labels = c("Normal","SN", "Gehan", "GEE", "NSNSM"), cex.axis = 5)
axis(2, cex.axis = 3)

######################################
######################c################################################################
## Assumption4. logT ~ Skewt(xi = 0, sigma = 1, lambda = -4, nu = 3)
######################################################################################
#### Generate Data #################################
###################################################

X = list()
logT = list()  # True log survival time
logY = list()  # Observed log survival time
delta = list() # delta = 1: observed, delta = 0: censored
Eps_true = list()

set.seed(9151321)
for(i in 1 : iter){

	# Covariate
		X_1 = rnorm(n_sample, 0, 1)
		X_2 = rbinom(n_sample, 1, 0.5)
		X[[i]] = cbind(X_1, X_2)
	
	#True logT
		Eps_true[[i]] = rst(n_sample, xi = 0, omega = 1, alpha = -4, nu = 3)
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
#L2 norm

L2_normal = NULL 
L2_skewnormal = NULL
L2_gehan = NULL
L2_gee = NULL
L2_NSNSM = NULL


for(i in 1 : iter){
	L2_normal = append(L2_normal, sqrt(sum((beta_normal[[i]] - beta_true)^2)))
	L2_skewnormal = append(L2_skewnormal, sqrt(sum((beta_skewnormal[[i]] - beta_true)^2)))
	L2_gehan = append(L2_gehan, sqrt(sum((beta_gehan[[i]] - beta_true)^2)))
	L2_gee = append(L2_gee, sqrt(sum((beta_gee[[i]] - beta_true)^2)))
	L2_NSNSM = append(L2_NSNSM, sqrt(sum((beta_NSNSM[[i]] - beta_true)^2)))
}


boxplot(L2_normal, L2_skewnormal, L2_gehan, L2_gee, L2_NSNSM, 
		yaxt = "n", col = c(rep("grey",4), rep("pink", 1)), ylim = c(0,2.5))
axis(1, at = 1:5, labels = c("Normal","SN", "Gehan", "GEE", "NSNSM"), cex.axis = 5)
axis(2, cex.axis = 3)
######################################

