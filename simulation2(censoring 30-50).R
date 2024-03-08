###########################################################################################
################## simulation 3 (Prediction) ###########################################
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
n_sample = 250
beta_true = c(2, 1, -1) 
tau = 4

test_sample = 1000


######################c################################################################
## Assumption1. logT ~ Skewt(xi = 0, sigma = 1, lambda = -15, nu = 3)
######################################################################################
#### Generate Data #################################
###################################################
X = list()
X_test = list()
logT = list()  # True log survival time
logT_test = list()
logY = list()  # Observed log survival time
logY_test = list()
delta = list() # delta = 1: observed, delta = 0: censored
delta_test = list()
Eps_true = list()
Eps_true_test = list()

set.seed(1234)
for(i in 1 : iter){

	# Covariate (train)
		X_1 = rnorm(n_sample, 0, 1)
		X_2 = rbinom(n_sample, 1, 0.5)
		X[[i]] = cbind(X_1, X_2)

	# Covariate (test)
		X_1 = rnorm(test_sample, 0, 1)
		X_2 = rbinom(test_sample, 1, 0.5)
		X_test[[i]] = cbind(X_1, X_2)

	
	#True logT (train)
		Eps_true[[i]] = rst(n_sample, xi = 0, omega = 1, alpha = -15, nu = 3)
		Eps_true[[i]] = (Eps_true[[i]] - mean(Eps_true[[i]]))/sd(Eps_true[[i]])

	   design.X = cbind(1, X[[i]])
		logT[[i]] = (design.X %*% beta_true) + Eps_true[[i]]

	#True logT (test)
		Eps_true_test[[i]] = rst(test_sample, xi = 0, omega = 1, alpha = -15, nu = 3)
		Eps_true_test[[i]] = (Eps_true_test[[i]] - mean(Eps_true_test[[i]]))/sd(Eps_true_test[[i]])

	   design.X_test = cbind(1, X_test[[i]])
		logT_test[[i]] = (design.X_test %*% beta_true) + Eps_true_test[[i]]

	#Censor: delta (train)
		logC <- runif(n_sample, 0, tau)	

		delta[[i]] = (logT[[i]] < logC) * 1	

	#Censor: delta (test)
		logC_test <- runif(test_sample, 0, tau)	

		delta_test[[i]] = (logT_test[[i]] < logC_test) * 1	

	# censored data (train)
		logY[[i]] = pmin(logT[[i]], logC)  #logT[[i]]

		if(delta[[i]][which.max(Eps_true[[i]]),] == 0){
			logY[[i]][which.max(Eps_true[[i]]),] = logT[[i]][which.max(Eps_true[[i]]),]
			delta[[i]][which.max(Eps_true[[i]]),] = 1
		}

	# censored data (test)
		logY_test[[i]] = pmin(logT_test[[i]], logC_test)  #logT[[i]]

		if(delta_test[[i]][which.max(Eps_true_test[[i]]),] == 0){
			logY_test[[i]][which.max(Eps_true_test[[i]]),] = logT_test[[i]][which.max(Eps_true_test[[i]]),]
			delta_test[[i]][which.max(Eps_true_test[[i]]),] = 1
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
gumbel = list()	   ; beta_gumbel = list()
skewnormal = list() ; beta_skewnormal = list()
NSNSM = list()      ;  beta_NSNSM = list()

RMSE_normal = NULL
RMSE_NSNSM = NULL
RMSE_SN = NULL

iter = 200
for(i in 1:iter){

	# train
	x = X[[i]]
	y = logY[[i]]
	d = delta[[i]]

	# test
	x_test = X_test[[i]]
	logt_test = logT_test[[i]]

	# training 
	normal[[i]] = survreg(Surv(time = exp(y), event = d) ~ x, dist = 'lognormal')
	skewnormal[[i]] = AFT_SN(x = x, y = y, delta = d)
	NSNSM[[i]] = AFT_NSNSM(x, y, d, simulation_index = i)

	beta_normal[[i]]  = normal[[i]]$coefficients %>% unname()
	beta_skewnormal[[i]]  = skewnormal[[i]][[1]]
	beta_NSNSM[[i]]  = NSNSM[[i]][[1]]

	#prediction
	RMSE_normal = append(RMSE_normal, sqrt(sum((logt_test - cbind(1, x_test) %*% beta_normal[[i]])^2)/test_sample))
	RMSE_NSNSM = append(RMSE_NSNSM, sqrt(sum((logt_test - cbind(1, x_test) %*% beta_NSNSM[[i]])^2)/test_sample))
	RMSE_SN = append(RMSE_SN, sqrt(sum((logt_test - cbind(1, x_test) %*% beta_skewnormal[[i]])^2)/test_sample))

	print(i)
}


### Gehan ###
gehan = list()
beta_gehan = list()

RMSE_gehan = NULL

for(i in 1 : iter){

	# train
	x = X[[i]]
	y = logY[[i]]
	d = delta[[i]]

	# test
	x_test = X_test[[i]]
	logt_test = logT_test[[i]]

	# training
	dt = as.data.frame(cbind(exp(y), x, d))
	colnames(dt) = c("Y", "X1", "X2", "delta")
	gehan[[i]] =  aftsrr(Surv(Y, delta) ~ X1 + X2, data = dt)
	beta_gehan[[i]] = gehan[[i]]$beta

	esp = y - x %*% as.matrix(gehan[[i]]$beta)
	
	KME <- survfit(Surv(esp, d)~ 1, conf.type="none")
	intercept = KME$time%*%diff(c(0,(1-KME$surv))) # Estimated intercept for Gehan

	beta_gehan[[i]] = c(intercept, beta_gehan[[i]])

	#prediction
	RMSE_gehan = append(RMSE_gehan, sqrt(sum((logt_test - cbind(1, x_test) %*% beta_gehan[[i]])^2)/test_sample)	)

	print(i)
}


## GEE ####
gee = list()
beta_gee = list()

RMSE_GEE = NULL

for(i in 1 : iter){

	# train
	x = X[[i]]
	y = logY[[i]]
	d = delta[[i]]

	# test
	x_test = X_test[[i]]
	logt_test = logT_test[[i]]

	# training
	dt = as.data.frame(cbind(exp(y), x, d))
	colnames(dt) = c("Y", "X1", "X2", "delta")
	gee[[i]] =  aftgee(Surv(Y, delta) ~ X1 + X2, data = dt)
	beta_gee[[i]] = gee[[i]]$coef.res

	#prediction
	RMSE_GEE = append(RMSE_GEE, sqrt(sum((logt_test - cbind(1, x_test) %*% beta_gee[[i]])^2)/test_sample)	)

	print(i)
}


######################################
#RMSE

boxplot(RMSE_normal, RMSE_SN, RMSE_gehan, RMSE_GEE, RMSE_NSNSM,  
		yaxt = "n", col = c(rep("grey",4), rep("pink", 1)), ylim = c(0.995,1.07))
axis(1, at = 1:5, labels = c("Normal", "SN", "Gehan", "GEE", "SSNSM"), cex.axis = 5)
axis(2, cex.axis = 3)


###########################################################
######################################################################################
## Assumption2. logT ~ Normal(mean = 0, sd = 1)
######################################################################################
#### Generate Data #################################
###################################################
X2 = list()
X2_test = list()
logT2 = list()  # True log survival time
logT2_test = list()
logY2 = list()  # Observed log survival time
logY2_test = list()
delta2 = list() # delta = 1: observed, delta = 0: censored
delta2_test = list()
Eps_true2 = list()
Eps_true2_test = list()

set.seed(1234)
for(i in 1 : iter){

	# Covariate (train)
		X_1 = rnorm(n_sample, 0, 1)
		X_2 = rbinom(n_sample, 1, 0.5)
		X2[[i]] = cbind(X_1, X_2)

	# Covariate (test)
		X_1 = rnorm(test_sample, 0, 1)
		X_2 = rbinom(test_sample, 1, 0.5)
		X2_test[[i]] = cbind(X_1, X_2)

	
	#True logT (train)
		Eps_true2[[i]] = rnorm(n_sample, mean = 0, sd = 1)
#		Eps_true2[[i]] = (Eps_true2[[i]] - mean(Eps_true2[[i]]))/sd(Eps_true2[[i]])

	   design.X2 = cbind(1, X2[[i]])
		logT2[[i]] = (design.X2 %*% beta_true) + Eps_true2[[i]]

	#True logT (test)
		Eps_true2_test[[i]] = rnorm(test_sample, mean = 0, sd = 1)
#		Eps_true2_test[[i]] = (Eps_true2_test[[i]] - mean(Eps_true2_test[[i]]))/sd(Eps_true2_test[[i]])

	   design.X2_test = cbind(1, X2_test[[i]])
		logT2_test[[i]] = (design.X2_test %*% beta_true) + Eps_true2_test[[i]]

	#Censor: delta (train)
		logC2 <- runif(n_sample, 0, tau)	

		delta2[[i]] = (logT2[[i]] < logC2) * 1	

	#Censor: delta (test)
		logC2_test <- runif(test_sample, 0, tau)	

		delta2_test[[i]] = (logT2_test[[i]] < logC2_test) * 1	

	# censored data (train)
		logY2[[i]] = pmin(logT2[[i]], logC2)  #logT[[i]]

		if(delta2[[i]][which.max(Eps_true2[[i]]),] == 0){
			logY2[[i]][which.max(Eps_true2[[i]]),] = logT2[[i]][which.max(Eps_true2[[i]]),]
			delta2[[i]][which.max(Eps_true2[[i]]),] = 1
		}

	# censored data (test)
		logY2_test[[i]] = pmin(logT2_test[[i]], logC2_test)  #logT[[i]]

		if(delta2_test[[i]][which.max(Eps_true2_test[[i]]),] == 0){
			logY2_test[[i]][which.max(Eps_true2_test[[i]]),] = logT2_test[[i]][which.max(Eps_true2_test[[i]]),]
			delta2_test[[i]][which.max(Eps_true2_test[[i]]),] = 1
		}
		
}

## censored rate:
censored_rate = NULL
for(i in 1 : iter){
	censored_rate[i] = mean(delta2[[i]] == 0)
}
summary(censored_rate)   


####################################
### Results #######################
###################################

normal2 = list()     ; beta_normal2 = list()
NSNSM2 = list()      ;  beta_NSNSM2 = list()
skewnormal2 = list()      ;  beta_skewnormal2 = list()

RMSE_normal2 = NULL
RMSE_NSNSM2 = NULL
RMSE_SN2 = NULL

iter = 200
for(i in 1:iter){

	# train
	x = X2[[i]]
	y = logY2[[i]]
	d = delta2[[i]]

	# test
	x_test = X2_test[[i]]
	logt_test = logT2_test[[i]]

	# training 
	normal2[[i]] = survreg(Surv(time = exp(y), event = d) ~ x, dist = 'lognormal')
	NSNSM2[[i]] = AFT_NSNSM(x, y, d, simulation_index = i)
	skewnormal2[[i]] = AFT_SN(x = x, y = y, delta = d)

	beta_normal2[[i]]  = normal2[[i]]$coefficients %>% unname()
	beta_NSNSM2[[i]]  = NSNSM2[[i]][[1]]
	beta_skewnormal2[[i]]  = skewnormal2[[i]][[1]]

	#prediction
	RMSE_normal2 = append(RMSE_normal2, sqrt(sum((logt_test - cbind(1, x_test) %*% beta_normal2[[i]])^2)/test_sample))
	RMSE_NSNSM2 = append(RMSE_NSNSM2, sqrt(sum((logt_test - cbind(1, x_test) %*% beta_NSNSM2[[i]])^2)/test_sample))
	RMSE_SN2 = append(RMSE_SN2, sqrt(sum((logt_test - cbind(1, x_test) %*% beta_skewnormal2[[i]])^2)/test_sample))

	print(i)
}



### Gehan ###
gehan2 = list()
beta_gehan2 = list()

RMSE_gehan2 = NULL

for(i in 1 : iter){

	# train
	x = X2[[i]]
	y = logY2[[i]]
	d = delta2[[i]]

	# test
	x_test = X2_test[[i]]
	logt_test = logT2_test[[i]]

	# training
	dt = as.data.frame(cbind(exp(y), x, d))
	colnames(dt) = c("Y", "X1", "X2", "delta")
	gehan2[[i]] =  aftsrr(Surv(Y, delta) ~ X1 + X2, data = dt)
	beta_gehan2[[i]] = gehan2[[i]]$beta

	esp = y - x %*% as.matrix(gehan2[[i]]$beta)
	
	KME <- survfit(Surv(esp, d)~ 1, conf.type="none")
	intercept = KME$time%*%diff(c(0,(1-KME$surv))) # Estimated intercept for Gehan

	beta_gehan2[[i]] = c(intercept, beta_gehan2[[i]])

	#prediction
	RMSE_gehan2 = append(RMSE_gehan2, sqrt(sum((logt_test - cbind(1, x_test) %*% beta_gehan2[[i]])^2)/test_sample)	)

	print(i)
}


## GEE ####
gee2 = list()
beta_gee2 = list()

RMSE_GEE2 = NULL

for(i in 1 : iter){

	# train
	x = X2[[i]]
	y = logY2[[i]]
	d = delta2[[i]]

	# test
	x_test = X2_test[[i]]
	logt_test = logT2_test[[i]]

	# training
	dt = as.data.frame(cbind(exp(y), x, d))
	colnames(dt) = c("Y", "X1", "X2", "delta")
	gee2[[i]] =  aftgee(Surv(Y, delta) ~ X1 + X2, data = dt)
	beta_gee2[[i]] = gee2[[i]]$coef.res

	#prediction
	RMSE_GEE2 = append(RMSE_GEE2, sqrt(sum((logt_test - cbind(1, x_test) %*% beta_gee2[[i]])^2)/test_sample)	)

	print(i)
}

######################################
#RMSE

boxplot(RMSE_normal2, RMSE_SN2, RMSE_gehan2, RMSE_GEE2, RMSE_NSNSM2,  
		yaxt = "n", col = c(rep("grey",4), rep("pink", 1)))
axis(1, at = 1:5, labels = c("Normal","SN","Gehan", "GEE", "SSNSM"), cex.axis = 5)
axis(2, cex.axis = 3)

###########################################################
######################################################################################
## Assumption3. logT ~ t(df = 3)
######################################################################################
#### Generate Data #################################
###################################################
X3 = list()
X3_test = list()
logT3 = list()  # True log survival time
logT3_test = list()
logY3 = list()  # Observed log survival time
logY3_test = list()
delta3 = list() # delta = 1: observed, delta = 0: censored
delta3_test = list()
Eps_true3 = list()
Eps_true3_test = list()

set.seed(1334)
for(i in 1 : iter){

	# Covariate (train)
		X_1 = rnorm(n_sample, 0, 1)
		X_2 = rbinom(n_sample, 1, 0.5)
		X3[[i]] = cbind(X_1, X_2)

	# Covariate (test)
		X_1 = rnorm(test_sample, 0, 1)
		X_2 = rbinom(test_sample, 1, 0.5)
		X3_test[[i]] = cbind(X_1, X_2)

	
	#True logT (train)
		Eps_true3[[i]] = rt(n_sample, df = 3)
		Eps_true3[[i]] = (Eps_true3[[i]] - mean(Eps_true3[[i]]))/sd(Eps_true3[[i]])

	   design.X3 = cbind(1, X3[[i]])
		logT3[[i]] = (design.X3 %*% beta_true) + Eps_true3[[i]]

	#True logT (test)
		Eps_true3_test[[i]] = rt(test_sample, df =3)
		Eps_true3_test[[i]] = (Eps_true3_test[[i]] - mean(Eps_true3_test[[i]]))/sd(Eps_true3_test[[i]])

	   design.X3_test = cbind(1, X3_test[[i]])
		logT3_test[[i]] = (design.X3_test %*% beta_true) + Eps_true3_test[[i]]

	#Censor: delta (train)
		logC3 <- runif(n_sample, 0, tau)	

		delta3[[i]] = (logT3[[i]] < logC3) * 1	

	#Censor: delta (test)
		logC3_test <- runif(test_sample, 0, tau)	

		delta3_test[[i]] = (logT3_test[[i]] < logC3_test) * 1	

	# censored data (train)
		logY3[[i]] = pmin(logT3[[i]], logC3)  #logT[[i]]

		if(delta3[[i]][which.max(Eps_true3[[i]]),] == 0){
			logY3[[i]][which.max(Eps_true3[[i]]),] = logT3[[i]][which.max(Eps_true3[[i]]),]
			delta3[[i]][which.max(Eps_true3[[i]]),] = 1
		}

	# censored data (test)
		logY3_test[[i]] = pmin(logT3_test[[i]], logC3_test)  #logT[[i]]

		if(delta3_test[[i]][which.max(Eps_true3_test[[i]]),] == 0){
			logY3_test[[i]][which.max(Eps_true3_test[[i]]),] = logT3_test[[i]][which.max(Eps_true3_test[[i]]),]
			delta3_test[[i]][which.max(Eps_true3_test[[i]]),] = 1
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
NSNSM3 = list()      ;  beta_NSNSM3 = list()
skewnormal3 = list()      ;  beta_skewnormal3 = list()

RMSE_normal3 = NULL
RMSE_NSNSM3 = NULL
RMSE_SN3 = NULL

iter = 200
for(i in 1:iter){

	# train
	x = X3[[i]]
	y = logY3[[i]]
	d = delta3[[i]]

	# test
	x_test = X3_test[[i]]
	logt_test = logT3_test[[i]]

	# training 
	normal3[[i]] = survreg(Surv(time = exp(y), event = d) ~ x, dist = 'lognormal')
	NSNSM3[[i]] = AFT_NSNSM(x, y, d, simulation_index = i)
	skewnormal3[[i]] = AFT_SN(x = x, y = y, delta = d)

	beta_normal3[[i]]  = normal3[[i]]$coefficients %>% unname()
	beta_NSNSM3[[i]]  = NSNSM3[[i]][[1]]
	beta_skewnormal3[[i]]  = skewnormal3[[i]][[1]]

	#prediction
	RMSE_normal3 = append(RMSE_normal3, sqrt(sum((logt_test - cbind(1, x_test) %*% beta_normal3[[i]])^2)/test_sample))
	RMSE_NSNSM3 = append(RMSE_NSNSM3, sqrt(sum((logt_test - cbind(1, x_test) %*% beta_NSNSM3[[i]])^2)/test_sample))
	RMSE_SN3 = append(RMSE_SN3, sqrt(sum((logt_test - cbind(1, x_test) %*% beta_skewnormal3[[i]])^2)/test_sample))

	print(i)
}


### Gehan ###
gehan3 = list()
beta_gehan3 = list()

RMSE_gehan3 = NULL

for(i in 1 : iter){

	# train
	x = X3[[i]]
	y = logY3[[i]]
	d = delta3[[i]]

	# test
	x_test = X3_test[[i]]
	logt_test = logT3_test[[i]]

	# training
	dt = as.data.frame(cbind(exp(y), x, d))
	colnames(dt) = c("Y", "X1", "X2", "delta")
	gehan3[[i]] =  aftsrr(Surv(Y, delta) ~ X1 + X2, data = dt)
	beta_gehan3[[i]] = gehan3[[i]]$beta

	esp = y - x %*% as.matrix(gehan3[[i]]$beta)
	
	KME <- survfit(Surv(esp, d)~ 1, conf.type="none")
	intercept = KME$time%*%diff(c(0,(1-KME$surv))) # Estimated intercept for Gehan

	beta_gehan3[[i]] = c(intercept, beta_gehan3[[i]])

	#prediction
	RMSE_gehan3 = append(RMSE_gehan3, sqrt(sum((logt_test - cbind(1, x_test) %*% beta_gehan3[[i]])^2)/test_sample)	)

	print(i)
}


## GEE ####
gee3 = list()
beta_gee3 = list()

RMSE_GEE3 = NULL

for(i in 1 : iter){

	# train
	x = X3[[i]]
	y = logY3[[i]]
	d = delta3[[i]]

	# test
	x_test = X3_test[[i]]
	logt_test = logT3_test[[i]]

	# training
	dt = as.data.frame(cbind(exp(y), x, d))
	colnames(dt) = c("Y", "X1", "X2", "delta")
	gee3[[i]] =  aftgee(Surv(Y, delta) ~ X1 + X2, data = dt)
	beta_gee3[[i]] = gee3[[i]]$coef.res

	#prediction
	RMSE_GEE3 = append(RMSE_GEE3, sqrt(sum((logt_test - cbind(1, x_test) %*% beta_gee3[[i]])^2)/test_sample)	)

	print(i)
}

######################################
#RMSE

boxplot(RMSE_normal3, RMSE_SN3, RMSE_gehan3, RMSE_GEE3, RMSE_NSNSM3,  
		yaxt = "n", ylim = c(0.99, 1.04) ,col = c(rep("grey",4), rep("pink", 1)))
axis(1, at = 1:5, labels = c("Normal","SN","Gehan", "GEE", "SSNSM"), cex.axis = 5)
axis(2, cex.axis = 3)


######################################################################################
## Assumption4. logT ~ Gumbel(0,5)  rgumbel(n_sample, loc = 0, scale = 1)
######################################################################################
#### Generate Data #################################
###################################################
X4 = list()
X4_test = list()
logT4 = list()  # True log survival time
logT4_test = list()
logY4 = list()  # Observed log survival time
logY4_test = list()
delta4 = list() # delta = 1: observed, delta = 0: censored
delta4_test = list()
Eps_true4 = list()
Eps_true4_test = list()

set.seed(1444)
for(i in 1 : iter){
	library(evd)
	# Covariate (train)
		X_1 = rnorm(n_sample, 0, 1)
		X_2 = rbinom(n_sample, 1, 0.5)
		X4[[i]] = cbind(X_1, X_2)

	# Covariate (test)
		X_1 = rnorm(test_sample, 0, 1)
		X_2 = rbinom(test_sample, 1, 0.5)
		X4_test[[i]] = cbind(X_1, X_2)

	
	#True logT (train)
		Eps_true4[[i]] = rgumbel(n_sample, loc = 0, scale = 5)
		Eps_true4[[i]] = (Eps_true4[[i]] - mean(Eps_true4[[i]]))/sd(Eps_true4[[i]])

	   design.X4 = cbind(1, X4[[i]])
		logT4[[i]] = (design.X4 %*% beta_true) + Eps_true4[[i]]

	#True logT (test)
		Eps_true4_test[[i]] = rgumbel(test_sample, loc = 0, scale = 5)
		Eps_true4_test[[i]] = (Eps_true4_test[[i]] - mean(Eps_true4_test[[i]]))/sd(Eps_true4_test[[i]])

	   design.X4_test = cbind(1, X4_test[[i]])
		logT4_test[[i]] = (design.X4_test %*% beta_true) + Eps_true4_test[[i]]

	#Censor: delta (train)
		logC4 <- runif(n_sample, 0, tau)	

		delta4[[i]] = (logT4[[i]] < logC4) * 1	

	#Censor: delta (test)
		logC4_test <- runif(test_sample, 0, tau)	

		delta4_test[[i]] = (logT4_test[[i]] < logC4_test) * 1	

	# censored data (train)
		logY4[[i]] = pmin(logT4[[i]], logC4)  #logT[[i]]

		if(delta4[[i]][which.max(Eps_true4[[i]]),] == 0){
			logY4[[i]][which.max(Eps_true4[[i]]),] = logT4[[i]][which.max(Eps_true4[[i]]),]
			delta4[[i]][which.max(Eps_true4[[i]]),] = 1
		}

	# censored data (test)
		logY4_test[[i]] = pmin(logT4_test[[i]], logC4_test)  #logT[[i]]

		if(delta4_test[[i]][which.max(Eps_true4_test[[i]]),] == 0){
			logY4_test[[i]][which.max(Eps_true4_test[[i]]),] = logT4_test[[i]][which.max(Eps_true4_test[[i]]),]
			delta4_test[[i]][which.max(Eps_true4_test[[i]]),] = 1
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
NSNSM4 = list()      ;  beta_NSNSM4 = list()
skewnormal4 = list()      ;  beta_skewnormal4 = list()

RMSE_normal4 = NULL
RMSE_NSNSM4 = NULL
RMSE_SN4 = NULL

iter = 200
for(i in 1:iter){

	# train
	x = X4[[i]]
	y = logY4[[i]]
	d = delta4[[i]]

	# test
	x_test = X4_test[[i]]
	logt_test = logT4_test[[i]]

	# training 
	normal4[[i]] = survreg(Surv(time = exp(y), event = d) ~ x, dist = 'lognormal')
	NSNSM4[[i]] = AFT_NSNSM(x, y, d, simulation_index = i)
	skewnormal4[[i]] = AFT_SN(x = x, y = y, delta = d)

	beta_normal4[[i]]  = normal4[[i]]$coefficients %>% unname()
	beta_NSNSM4[[i]]  = NSNSM4[[i]][[1]]
	beta_skewnormal4[[i]]  = skewnormal4[[i]][[1]]

	#prediction
	RMSE_normal4 = append(RMSE_normal4, sqrt(sum((logt_test - cbind(1, x_test) %*% beta_normal4[[i]])^2)/test_sample))
	RMSE_NSNSM4 = append(RMSE_NSNSM4, sqrt(sum((logt_test - cbind(1, x_test) %*% beta_NSNSM4[[i]])^2)/test_sample))
	RMSE_SN4 = append(RMSE_SN4, sqrt(sum((logt_test - cbind(1, x_test) %*% beta_skewnormal4[[i]])^2)/test_sample))

	print(i)
}


### Gehan ###
gehan4 = list()
beta_gehan4 = list()

RMSE_gehan4 = NULL

for(i in 1 : iter){

	# train
	x = X4[[i]]
	y = logY4[[i]]
	d = delta4[[i]]

	# test
	x_test = X4_test[[i]]
	logt_test = logT4_test[[i]]

	# training
	dt = as.data.frame(cbind(exp(y), x, d))
	colnames(dt) = c("Y", "X1", "X2", "delta")
	gehan4[[i]] =  aftsrr(Surv(Y, delta) ~ X1 + X2, data = dt)
	beta_gehan4[[i]] = gehan4[[i]]$beta

	esp = y - x %*% as.matrix(gehan4[[i]]$beta)
	
	KME <- survfit(Surv(esp, d)~ 1, conf.type="none")
	intercept = KME$time%*%diff(c(0,(1-KME$surv))) # Estimated intercept for Gehan

	beta_gehan4[[i]] = c(intercept, beta_gehan4[[i]])

	#prediction
	RMSE_gehan4 = append(RMSE_gehan4, sqrt(sum((logt_test - cbind(1, x_test) %*% beta_gehan4[[i]])^2)/test_sample)	)

	print(i)
}


## GEE ####
gee4 = list()
beta_gee4 = list()

RMSE_GEE4 = NULL

for(i in 1 : iter){

	# train
	x = X4[[i]]
	y = logY4[[i]]
	d = delta4[[i]]

	# test
	x_test = X4_test[[i]]
	logt_test = logT4_test[[i]]

	# training
	dt = as.data.frame(cbind(exp(y), x, d))
	colnames(dt) = c("Y", "X1", "X2", "delta")
	gee4[[i]] =  aftgee(Surv(Y, delta) ~ X1 + X2, data = dt)
	beta_gee4[[i]] = gee4[[i]]$coef.res

	#prediction
	RMSE_GEE4 = append(RMSE_GEE4, sqrt(sum((logt_test - cbind(1, x_test) %*% beta_gee4[[i]])^2)/test_sample)	)

	print(i)
}


######################################
#RMSE

boxplot(RMSE_normal4, RMSE_SN4, RMSE_gehan4, RMSE_GEE4, RMSE_NSNSM4,  
		yaxt = "n", ylim = c(0.995,1.03), col = c(rep("grey",4), rep("pink", 1)))
axis(1, at = 1:5, labels = c("Normal", "SN", "Gehan", "GEE", "SSNSM"), cex.axis = 5)
axis(2, cex.axis = 3)

##############################################################