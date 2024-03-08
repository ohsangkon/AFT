############################################################
library(tidyverse)
library(sn)
library(nnls)
library(latex2exp)
library(e1071)
library(survival)
library(aftgee)
library(SurvMetrics)
library(stats) # ecdf

options(scipen = 100)

############################################################
# EDA
colnames(lung)
str(lung)

sum(which(apply(is.na(lung), 1, sum) != 0))

### complete data
lung_comp = na.omit(lung)
lung_comp$status = (lung_comp$status - 1)
lung_comp = lung_comp[,-1]

lung_comp$ph.karno = log(lung_comp$ph.karno)
lung_comp$pat.karno = log(lung_comp$pat.karno)
lung_comp$meal.cal = log(lung_comp$meal.cal)

#lung_comp$sex = factor(lung_comp$sex)

str(lung_comp)
summary(lung_comp)


y = as.matrix(lung_comp[,1])
x = as.matrix(lung_comp[,-c(1,2)])
d = lung_comp[,2]
prop.table(table(d)) # 28.1% censoring

t = lung_comp[lung_comp$status == 1,1]
par(mar=c(5,5,5,5))
hist(log(t), freq = F, cex.axis = 1.5, cex.lab = 2, xlab = "log(T)", main = "")


############################################################
## Bootstrap dataset

B = 500 # Bootstrap size

lung_comp_boots = list()

set.seed(7821)
boots_num = lapply(1:B, FUN = function(i){ 
						sample(1:nrow(lung_comp), 
						size = nrow(lung_comp), replace = TRUE)}
					)

lung_comp_boots = lapply(1:B, FUN = function (i){
						lung_comp[boots_num[[i]],]
					})


y2 = list()
for(i in 1 : B){y2[[i]] = as.matrix(lung_comp_boots[[i]]$time)}


x2 = list()
for(i in 1 : B){
	attach(lung_comp_boots[[i]])
	x2[[i]] = cbind(age, sex, ph.ecog, ph.karno, pat.karno, 
										    meal.cal, wt.loss)

	detach(lung_comp_boots[[i]])
}


d2 = list()
for(i in 1 : B){
	attach(lung_comp_boots[[i]])
	d2[[i]] = lung_comp_boots[[i]]$status==1
	detach(lung_comp_boots[[i]])
}
d2

##########################
# normal
normal = survreg(data = lung_comp, Surv(time = time, event = status) ~age + sex + ph.ecog + ph.karno + pat.karno + 
																								meal.cal + wt.loss, dist = 'lognormal')
normal$coefficients


#bootstrap
PAFT_boots = lapply(1:B, FUN = function(i){
							survreg(Surv(time = time, event = status) ~age + sex + ph.ecog + ph.karno + pat.karno + 
																								meal.cal + wt.loss
				, data=lung_comp_boots[[i]] ,dist="lognormal")
							}
					)

PAFT_boots_coeff = sapply(1:B, FUN = function(i){
							PAFT_boots[[i]]$coefficient
							}
					) 

normal_sd = sqrt(apply((PAFT_boots_coeff - apply(PAFT_boots_coeff,1,mean))^2,1,sum)/ (B-1))


###########################
# gee
gee =  aftgee(Surv(time, status) ~ age + sex + ph.ecog + ph.karno + pat.karno + 
										    meal.cal + wt.loss, data = lung_comp)
gee$coef.res


###
# Bootstrap

SAFT_G_boots = lapply(1:B, FUN = function(i){
				aftgee(Surv(time, status) ~ age + sex + ph.ecog + ph.karno + pat.karno + 
										    meal.cal + wt.loss, 
					data=lung_comp_boots[[i]],B=0)
					}
				)

SAFT_G_boots_coeff = sapply(1:B, FUN = function(i){
							SAFT_G_boots[[i]]$coef.res							}
					) 

gee_sd = sqrt(apply((SAFT_G_boots_coeff - apply(SAFT_G_boots_coeff,1,mean))^2,1,sum)/ (B-1))

############################
# gehan
gehan =  aftsrr(Surv(time, status) ~ age + sex + ph.ecog + ph.karno + pat.karno + 
										    meal.cal + wt.loss, data = lung_comp)
gehan$beta

esp = log(y) - x %*% as.matrix(gehan$beta)
	
KME <- survfit(Surv(esp, d)~ 1, conf.type="none")
intercept = KME$time%*%diff(c(0,(1-KME$surv))) # Estimated intercept for Gehan

beta_gehan = c(intercept, gehan$beta)
beta_gehan


###
# Bootstrap
gehan_func = function(lung_comp_boots){

		SAFT_R = aftsrr(Surv(time, status) ~ age + sex + ph.ecog + ph.karno + pat.karno + 
										    meal.cal + wt.loss, 
						data=lung_comp_boots, B=0,se = "ISMB") #Gehan

		######## Calculation of the intercept
		attach(lung_comp_boots)
		data.x=as.matrix(cbind(age, sex, ph.ecog, ph.karno, pat.karno, 
										    meal.cal, wt.loss))
		detach(lung_comp_boots)

		mydata1 = log(lung_comp_boots$time) - data.x%*%as.matrix(SAFT_R$beta)
		ind=ifelse(lung_comp_boots$status==1,1,0)
		mydata.surv <- survfit(Surv(mydata1, ind)~ 1, conf.type="none")
		est.int = mydata.surv$time%*%diff(c(0,(1-mydata.surv$surv))) # Estimated intercept for Gehan

		SAFT_R_beta = c(est.int,SAFT_R$beta)
		names(SAFT_R_beta) = c("intercept", "age", "sex", "ph.ecog", "ph.karno",
										"pat.karno", "meal.cal", "wt.loss")
			
		return(SAFT_R_beta)

}

SAFT_R_boots_coeff = sapply(1:B, FUN = function(i){
							gehan_func(lung_comp_boots[[i]])
							}
					)


gehan_sd = sqrt(apply((SAFT_R_boots_coeff - apply(SAFT_R_boots_coeff,1,mean))^2,1,sum)/ (B-1))

######################################################
# SN
SN = AFT_SN(x = x, y = log(y), delta = d)

names(SN$beta_hat) = c("intercept", "age", "sex", "ph.ecog", "ph.karno",
										"pat.karno", "meal.cal", "wt.loss")

SN$beta_hat


##bootstrap
SN_boots = lapply(1:B, FUN = function(i){
							AFT_SN(x = x2[[i]], y = log(y2[[i]]), delta = d2[[i]])
							}
					)

SN_boots_coeff = sapply(1:B, FUN = function(i){
							SN_boots[[i]]$beta_hat
							}
					) 

SN_sd = sqrt(apply((SN_boots_coeff - apply(SN_boots_coeff,1,mean))^2,1,sum)/ (B-1))
SN_sd

##############################
# NSNSM
NSNSM = AFT_NSNSM(x = x, y = log(y), delta = d)
NSNSM$beta_hat

# Survival plot (Sex) 

X = lung_comp %>% select(-(time:status))
fixed = apply(X, 2, median)

for(i in 1:2){
	sex = fixed
	sex[2] = i

	Survival_prob = NULL 
	points = sort(lung_comp$time[which(d == 1)])

	Survival_prob = survival_prob_NSNSM(points, NSNSM, sex)

	Survival = cbind((points), Survival_prob)

	if(i == 1){
		par(mar = c(5,5,5,5))
		plot(Survival, type = "s", col = "blue", lwd = 2, cex.lab = 2,
			xlab = "Failure Time (Days)", ylab = "Conditional Survival Probability")
		legend("topright", legend = c("Male", "Female"), lty = 1, lwd = 2,
			col = c("blue", "red"))

	}else{
		lines(Survival, type = "s", col = "red", lwd = 2)
	}
}

# Survival plot (ECOG) 

fixed = apply(X, 2, median)

for(i in 1:4){
	ECOG = fixed
	ECOG[3] = i

	Survival_prob = NULL 
	points = sort(lung_comp$time[which(d == 1)])

	Survival_prob = survival_prob_NSNSM(points, NSNSM, ECOG)

	Survival = cbind((points), Survival_prob)

	if(i == 1){
		par(mar = c(5,5,5,5))
		plot(Survival, type = "s", col = "black", lwd = 2, cex.lab = 2,
			xlab = "Failure Time (Days)", ylab = "Conditional Survival Probability")
		legend("topright", legend = c("ECOG = 1", "ECOG = 2", "ECOG = 3", "ECOG = 4"), lty = 1, lwd = 2,
			col = c("black", "red", "blue", "chartreuse2"))

	}else if(i == 2){
		lines(Survival, type = "s", col = "red", lwd = 2)
	}else if(i == 3){
		lines(Survival, type = "s", col = "blue", lwd = 2)
	}else{
		lines(Survival, type = "s", col = "chartreuse2", lwd = 2)
	}
}


##bootstrap
NSNSM_boots = lapply(1:B, FUN = function(i){
							AFT_NSNSM(x = x2[[i]], y = log(y2[[i]]), delta = d2[[i]], simulation_index = i)
							}
					)

NSNSM_boots_coeff = sapply(1:B, FUN = function(i){
							NSNSM_boots[[i]]$beta_hat
							}
					) 

NSNSM_sd = sqrt(apply((NSNSM_boots_coeff - apply(NSNSM_boots_coeff,1,mean))^2,1,sum)/ (B-1))
NSNSM_sd


######################################################################
## Prediction (10-fold cross-validation) 

library(SurvMetrics)
library(caret)  
library(survival)  
library(pec)

######################################################################

index_data = list()
train_data = list()
test_data = list()

train_x = list()
train_y = list()
test_x = list()
d = list()

set.seed(530)
index_data = createFolds(1:nrow(lung_comp), 10)

for(i in 1 : 10){
 	 train_data[[i]] = lung_comp[-index_data[[i]],]
 	 test_data[[i]] = lung_comp[index_data[[i]],]

	attach(train_data[[i]])
	train_y[[i]] = as.matrix(train_data[[i]]$time)
	train_x[[i]] = cbind(age, sex, ph.ecog, ph.karno, pat.karno, 
										    meal.cal, wt.loss)
	d[[i]] = train_data[[i]]$status==1

	detach(train_data[[i]])

	attach(test_data[[i]])
	test_x[[i]] = cbind(age, sex, ph.ecog, ph.karno, pat.karno, 
										    meal.cal, wt.loss)

	detach(test_data[[i]])

}


#####################################################################
## Parametric MLE under normality
train_PAFT = list()
pred_PAFT = list()


for(i in 1 : 10){
	train_PAFT[[i]] = survreg(Surv(time,status==1)~ age + sex + ph.ecog + ph.karno + pat.karno + 
										    meal.cal + wt.loss
					, data=train_data[[i]] ,dist="lognormal")

}


#############################################################
# Fitting a smoothed Gehan estimator: using time in years

train_SAFT_R = list()
SAFT_R_beta = list()



for(i in 1 : 10){
	train_SAFT_R[[i]] = aftsrr(Surv(time,status==1) ~ age + sex + ph.ecog + ph.karno + pat.karno + 
										    meal.cal + wt.loss, 
						data=train_data[[i]], B=0,se = "ISMB") #Gehan

	attach(train_data[[i]])

	data.x=as.matrix(cbind(age, sex, ph.ecog, ph.karno, pat.karno, 
										    meal.cal, wt.loss))
	mydata1 = log(train_data[[i]]$time) - data.x%*%as.matrix(train_SAFT_R[[i]]$beta)
	ind=ifelse(train_data[[i]]$status==1,1,0)
	mydata.surv <- survfit(Surv(mydata1, ind)~ 1, conf.type="none")
	est.int = mydata.surv$time%*%diff(c(0,(1-mydata.surv$surv))) # Estimated intercept for Gehan

	SAFT_R_beta[[i]] = c(est.int,	train_SAFT_R[[i]]$beta)
	names(SAFT_R_beta[[i]]) = c("intercept", "age", "sex", "ph.ecog", "ph.karno", "pat.karno", 
										    "meal.cal", "wt.loss")

	detach(train_data[[i]])

}



####################################################################
# Least square estimator using GEE

train_SAFT_G = list()
pred_SAFT_G = list()


for(i in 1 : 10){
	train_SAFT_G[[i]] = aftgee(Surv(time, status==1)~ age + sex + ph.ecog + ph.karno + pat.karno + 
										    meal.cal + wt.loss, 
					data=train_data[[i]],B=0)
}

######################################################
### SN  ####

train_SN = list()

for(i in 1 : 10){
	attach(train_data[[i]])
	y = as.matrix(train_data[[i]]$time)
	x = cbind(age, sex, ph.ecog, ph.karno, pat.karno, 
										    meal.cal, wt.loss)
	d = train_data[[i]]$status==1

	detach(train_data[[i]])

	train_SN[[i]] = AFT_SN(x = x, y = log(y), delta = d)

}

######################################################
### Proposed method  ####

train_NSNSM = list()

for(i in 1 : 10){
	attach(train_data[[i]])
	y = as.matrix(train_data[[i]]$time)
	x = cbind(age, sex, ph.ecog, ph.karno, pat.karno, 
										    meal.cal, wt.loss)
	d = train_data[[i]]$status==1

	detach(train_data[[i]])

	train_NSNSM[[i]] = AFT_NSNSM(x = x, y = log(y), delta = d, simulation_index = i)

}


#####################################################################
#### IBS 

##
IBS_normal = NULL
IBS_NSNSM = NULL
IBS_GEE = NULL
IBS_Gehan = NULL
IBS_SN = NULL

for(i in 1 : 10){
	Surv_obj = Surv(test_data[[i]]$time, test_data[[i]]$status==1)

	pred_normal = NULL
	pred_NSNSM = NULL
	pred_Gehan = NULL
	pred_GEE = NULL
	pred_SN = NULL

	IBSrange = seq(from = 1, to = 365.25 * 5, length = 100)
	for(j in 1 : length(IBSrange)){
		pred_normal= cbind(pred_normal, survival_prob_Normal(t_star = IBSrange[j], model = train_PAFT[[i]], covariate = test_x[[i]]))
		pred_NSNSM= cbind(pred_NSNSM, survival_prob_NSNSM(t_star = IBSrange[j], model = train_NSNSM[[i]], covariate = test_x[[i]]))
		pred_GEE= cbind(pred_GEE, survival_prob_GEE(t_star = IBSrange[j], model = train_SAFT_G[[i]], covariate = test_x[[i]], train_x = train_x[[i]], train_y = train_y[[i]], d=d[[i]]))
		pred_Gehan= cbind(pred_Gehan, survival_prob_Gehan(t_star = IBSrange[j], model_beta = SAFT_R_beta[[i]], covariate = test_x[[i]], train_x = train_x[[i]], train_y= train_y[[i]], d=d[[i]]))
		pred_SN= cbind(pred_SN, survival_prob_SN(t_star = IBSrange[j], model = train_SN[[i]], covariate = test_x[[i]]))

	}
	IBS_normal = append( IBS_normal, IBS(Surv_obj, pred_normal, IBSrange = IBSrange))
	IBS_NSNSM = append( IBS_NSNSM, IBS(Surv_obj, pred_NSNSM, IBSrange = IBSrange))
	IBS_Gehan = append( IBS_Gehan, IBS(Surv_obj, pred_Gehan, IBSrange = IBSrange))
	IBS_GEE = append( IBS_GEE, IBS(Surv_obj, pred_GEE, IBSrange = IBSrange))
	IBS_SN = append( IBS_SN, IBS(Surv_obj, pred_SN, IBSrange = IBSrange))

	print(i)
}

mean(IBS_normal) ; mean(IBS_SN) ; mean(IBS_Gehan) ; mean(IBS_GEE) ; mean(IBS_NSNSM) 

