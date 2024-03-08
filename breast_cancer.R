#####################################################################
library(survival)
library(aftgee)
#library(rms)
library(tidyverse)
library(sn)
library(nnls)
library(latex2exp)
library(e1071)
library(emplik)
library(SurvMetrics)
library(stats) # ecdf

options(scipen=999)

#####################################################################
# Breast cancer dataset
############################################################################
breast = read.csv(".../Breast_Cancer.csv")
str(breast)

mybreast = breast %>% filter(Race == "Black")
mybreast$Status = factor(mybreast$Status)
prop.table(table(mybreast$Status))  # censoring rate: 67%

mybreast$Race = factor(mybreast$Race)
mybreast$Marital.Status = factor(mybreast$Marital.Status)
mybreast$T.Stage = factor(mybreast$T.Stage)
mybreast$N.Stage = factor(mybreast$N.Stage)
mybreast$X6th.Stage  = factor(mybreast$X6th.Stage)
mybreast$differentiate  = factor(mybreast$differentiate)
mybreast$Grade  = factor(mybreast$Grade)
mybreast$A.Stage = factor(mybreast$A.Stage)
mybreast$Estrogen.Status = factor(mybreast$Estrogen.Status)
mybreast$Progesterone.Status = factor(mybreast$Progesterone.Status)

str(mybreast)
summary(mybreast)

prop.table(table(mybreast$Status))  # censoring rate: 74.9%


attach(mybreast)
y = as.matrix(mybreast$Survival.Months)
x = cbind(Age, Tumor.Size, 
				Regional.Node.Examined, Reginol.Node.Positive)
d = mybreast$Status=="Dead"
prop.table(table(d)) # 61.5% censoring

detach(mybreast)

##################################################
### Bootstrap
###################################################

B = 500 # Bootstrap size

mybreast_boots = list()

set.seed(7821)
boots_num = lapply(1:B, FUN = function(i){ 
						sample(1:nrow(mybreast), 
						size = nrow(mybreast), replace = TRUE)}
					)

mybreast_boots = lapply(1:B, FUN = function (i){
						mybreast[boots_num[[i]],]
					})


y2 = list()
for(i in 1 : B){y2[[i]] = as.matrix(mybreast_boots[[i]]$Survival.Months)}


x2 = list()
for(i in 1 : B){
	attach(mybreast_boots[[i]])
	x2[[i]] = cbind(Age, Tumor.Size, 
				Regional.Node.Examined, Reginol.Node.Positive)

	detach(mybreast_boots[[i]])
}


d2 = list()
for(i in 1 : B){
	attach(mybreast_boots[[i]])
	d2[[i]] = mybreast_boots[[i]]$Status=="Dead"
	detach(mybreast_boots[[i]])
}
d2


###############################################################################
## Parametric MLE under normality
PAFT=survreg(Surv(Survival.Months,Status=="Dead")~Age + Tumor.Size + 
				 Regional.Node.Examined + Reginol.Node.Positive 
				, data=mybreast ,dist="lognormal")
PAFT$coefficient


#bootstrap
PAFT_boots = lapply(1:B, FUN = function(i){
							survreg(Surv(Survival.Months,Status=="Dead")~Age + Tumor.Size + 
				 Regional.Node.Examined + Reginol.Node.Positive
				, data=mybreast_boots[[i]] ,dist="lognormal")
							}
					)

PAFT_boots_coeff = sapply(1:B, FUN = function(i){
							PAFT_boots[[i]]$coefficient
							}
					) 

sqrt(apply((PAFT_boots_coeff - apply(PAFT_boots_coeff,1,mean))^2,1,sum)/ (B-1))

#############################################################
# Fitting a smoothed Gehan estimator: using Survival.Months in years
SAFT_R = aftsrr(Surv(Survival.Months,Status=="Dead")~Age + Tumor.Size + 
				 Regional.Node.Examined + Reginol.Node.Positive 
				, data=mybreast , B=0,se = "ISMB") #Gehan

######## Calculation of the intercept
attach(mybreast)
data.x=as.matrix(cbind(Age, Tumor.Size, 
				Regional.Node.Examined, Reginol.Node.Positive))
mydata1 = log(mybreast$Survival.Months) - data.x%*%as.matrix(SAFT_R$beta)
ind=ifelse(mybreast$Status=="Dead",1,0)
mydata.surv <- survfit(Surv(mydata1, ind)~ 1, conf.type="none")
est.int = mydata.surv$time%*%diff(c(0,(1-mydata.surv$surv))) # Estimated intercept for Gehan

SAFT_R_beta = c(est.int,SAFT_R$beta)
names(SAFT_R_beta) = c("intercept", "Age", "Tumor.Size", 
				"Regional.Node.Examined", "Reginol.Node.Positive")
SAFT_R_beta


detach(mybreast)

###
# Bootstrap
gehan_func = function(mybreast_boots){

		SAFT_R = aftsrr(Surv(Survival.Months,Status=="Dead") ~ Age + Tumor.Size + 
				 Regional.Node.Examined + Reginol.Node.Positive, 
						data=mybreast_boots, B=0,se = "ISMB") #Gehan

		######## Calculation of the intercept
		#attach(mybreast_boots)
		data.x=as.matrix(cbind(mybreast_boots$Age, mybreast_boots$Tumor.Size, 
				mybreast_boots$Regional.Node.Examined, mybreast_boots$Reginol.Node.Positive))
		#detach(mybreast_boots)

		mydata1 = log(mybreast_boots$Survival.Months) - data.x%*%as.matrix(SAFT_R$beta)
		ind=ifelse(mybreast_boots$Status=="Dead",1,0)
		mydata.surv <- survfit(Surv(mydata1, ind)~ 1, conf.type="none")
		est.int = mydata.surv$time%*%diff(c(0,(1-mydata.surv$surv))) # Estimated intercept for Gehan

		SAFT_R_beta = c(est.int,SAFT_R$beta)
		names(SAFT_R_beta) = c("intercept", "Age", "Tumor.Size", 
				"Regional.Node.Examined", "Reginol.Node.Positive")
			
		return(SAFT_R_beta)

}

SAFT_R_boots_coeff = sapply(1:B, FUN = function(i){
							gehan_func(mybreast_boots[[i]])
							}
					)


sqrt(apply((SAFT_R_boots_coeff - apply(SAFT_R_boots_coeff,1,mean))^2,1,sum)/ (B-1))

####################################################################
# Least square estimator using GEE
SAFT_G=aftgee(Surv(Survival.Months,Status=="Dead")~Age + Tumor.Size + 
				 Regional.Node.Examined + Reginol.Node.Positive 
				, data=mybreast ,B=0)

SAFT_G$coef.res

###
# Bootstrap

SAFT_G_boots = lapply(1:B, FUN = function(i){
				aftgee(Surv(Survival.Months, Status=="Dead")~ Age + Tumor.Size + 
				 Regional.Node.Examined + Reginol.Node.Positive, 
					data=mybreast_boots[[i]],B=0)
					}
				)

SAFT_G_boots_coeff = sapply(1:B, FUN = function(i){
							SAFT_G_boots[[i]]$coef.res							}
					) 

sqrt(apply((SAFT_G_boots_coeff - apply(SAFT_G_boots_coeff,1,mean))^2,1,sum)/ (B-1))

######################################################
### Skew normal  ####
attach(mybreast)
y = as.matrix(mybreast$Survival.Months)
x = cbind(Age, Tumor.Size, 
				Regional.Node.Examined, Reginol.Node.Positive)
d = mybreast$Status=="Dead"
prop.table(table(d)) # 61.5% censoring

detach(mybreast)

# SN
SN = AFT_SN(x = x, y = log(y), delta = d)

names(SN$beta_hat) = c("intercept", "Age", "Tumor.Size", 
				"Regional.Node.Examined", "Reginol.Node.Positive")

SN$beta_hat


#bootstrap
SN_boots = lapply(1:B, FUN = function(i){
							AFT_SN(x = x2[[i]], y = log(y2[[i]]), delta = d2[[i]])
							}
					)

SN_boots_coeff = sapply(1:B, FUN = function(i){
							SN_boots[[i]]$beta_hat
							}
					) 

sqrt(apply((SN_boots_coeff - apply(SN_boots_coeff,1,mean))^2,1,sum)/ (B-1))


######################################################
### Proposed method  ####
attach(mybreast)
y = as.matrix(mybreast$Survival.Months)
x = cbind(Age, Tumor.Size, 
				Regional.Node.Examined, Reginol.Node.Positive)
d = mybreast$Status=="Dead"
prop.table(table(d)) # 61.5% censoring

detach(mybreast)

# NSNSM
NSNSM = AFT_NSNSM(x = x, y = log(y), delta = d)

names(NSNSM$beta_hat) = c("intercept", "Age", "Tumor.Size", 
				"Regional.Node.Examined", "Reginol.Node.Positive")


NSNSM$beta_hat


### Survival plot for each case
X = mybreast %>% select(Age, Tumor.Size, 
				Regional.Node.Examined, Reginol.Node.Positive)

fixed = apply(X, 2, median)

positive = c(1, 5, 10, 20)

for(i in positive){
	covariate = fixed
	covariate[4] = i

	Survival_prob = NULL 
	points = sort(mybreast$Survival.Months[which(d == 1)])

	Survival_prob = survival_prob_NSNSM(points, NSNSM, covariate)

	Survival = cbind((points), Survival_prob)

	if(i == 1){
		par(mar = c(5,5,5,5))
		plot(Survival, type = "s", col = "black", lwd = 2, cex.lab = 2, ylim = c(0,1),
			xlab = "Failure Time (in months)", ylab = "Conditional Survival Probability")
		legend("bottomleft", legend = c("Node Positive = 1", "Node Positive = 5", "Node Positive = 10", "Node Positive = 20"), lty = 1, lwd = 2,
			col = c("black", "red", "blue", "chartreuse2"))

	}else if(i == 5){
		lines(Survival, type = "s", col = "red", lwd = 2)
	}else if(i == 10){
		lines(Survival, type = "s", col = "blue", lwd = 2)
	}else{
		lines(Survival, type = "s", col = "chartreuse2", lwd = 2)
	}
}


#bootstrap
NSNSM_boots = lapply(1:B, FUN = function(i){
							AFT_NSNSM(x = x2[[i]], y = log(y2[[i]]), delta = d2[[i]], simulation_index = i)
							}
					)

NSNSM_boots_coeff = sapply(1:B, FUN = function(i){
							NSNSM_boots[[i]]$beta_hat
							}
					) 

sd = sqrt(apply((NSNSM_boots_coeff - apply(NSNSM_boots_coeff,1,mean))^2,1,sum)/ (B-1))


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
index_data = createFolds(1:nrow(mybreast), 10)

for(i in 1 : 10){
 	 train_data[[i]] = mybreast[-index_data[[i]],]
 	 test_data[[i]] = mybreast[index_data[[i]],]

	attach(train_data[[i]])
	test_data[[i]]$Status = test_data[[i]]$Status 

	train_y[[i]] = as.matrix(train_data[[i]]$Survival.Months)
	train_x[[i]] = cbind(Age, Tumor.Size, 
				Regional.Node.Examined, Reginol.Node.Positive)
	d[[i]] = train_data[[i]]$Status=="Dead"

	detach(train_data[[i]])

	attach(test_data[[i]])
	test_x[[i]] = cbind(Age, Tumor.Size, 
				Regional.Node.Examined, Reginol.Node.Positive)

	detach(test_data[[i]])

}

#####################################################################
## Parametric MLE under normality
train_PAFT = list()

for(i in 1 : 10){
	train_PAFT[[i]] = survreg(Surv(Survival.Months, Status=="Dead")~ Age + Tumor.Size+ 
				Regional.Node.Examined+ Reginol.Node.Positive
					, data=train_data[[i]] ,dist="lognormal")

}


#############################################################
# Fitting a smoothed Gehan estimator: using Survival.Months in years

train_SAFT_R = list()
pred_SAFT_R = list()

SAFT_R_beta = list()


for(i in 1 : 10){
	train_SAFT_R[[i]] = aftsrr(Surv(Survival.Months, Status=="Dead")~ Age + Tumor.Size+ 
				Regional.Node.Examined+ Reginol.Node.Positive, 
						data=train_data[[i]], B=0,se = "ISMB") #Gehan

	attach(train_data[[i]])

	data.x=as.matrix(cbind(Age, Tumor.Size, 
				Regional.Node.Examined, Reginol.Node.Positive))
	mydata1 = log(train_data[[i]]$Survival.Months) - data.x%*%as.matrix(train_SAFT_R[[i]]$beta)
	ind=ifelse(train_data[[i]]$Status=="Dead",1,0)
	mydata.surv <- survfit(Surv(mydata1, ind)~ 1, conf.type="none")
	est.int = mydata.surv$time%*%diff(c(0,(1-mydata.surv$surv))) # Estimated intercept for Gehan

	SAFT_R_beta[[i]] = c(est.int,	train_SAFT_R[[i]]$beta)
	names(SAFT_R_beta[[i]]) = c("intercept", "Age", "Tumor.Size", 
				"Regional.Node.Examined", "Reginol.Node.Positive")

	detach(train_data[[i]])

}

####################################################################
# Least square estimator using GEE

train_SAFT_G = list()

for(i in 1 : 10){
	train_SAFT_G[[i]] = aftgee(Surv(Survival.Months, Status=="Dead")~ Age + Tumor.Size+ 
				Regional.Node.Examined+ Reginol.Node.Positive, 
						data=train_data[[i]],B=0)

}


######################################################
### SN  ####

train_SN = list()
pred_SN = list()


for(i in 1 : 10){
	attach(train_data[[i]])
	y = as.matrix(train_data[[i]]$Survival.Months)
	x = cbind(Age, Tumor.Size, 
				Regional.Node.Examined, Reginol.Node.Positive)
	d = train_data[[i]]$Status=="Dead"

	detach(train_data[[i]])

	train_SN[[i]] = AFT_SN(x = x, y = log(y), delta = d)

}


######################################################
### Proposed method  ####

train_NSNSM = list()

for(i in 1 : 10){
	attach(train_data[[i]])
	y = as.matrix(train_data[[i]]$Survival.Months)
	x = cbind(Age, Tumor.Size, 
				Regional.Node.Examined, Reginol.Node.Positive)
	d = train_data[[i]]$Status=="Dead"

	detach(train_data[[i]])

	train_NSNSM[[i]] = AFT_NSNSM(x = x, y = log(y), delta = d)

	print(i)

}


####################################
# IBS
####

IBS_normal = NULL
IBS_NSNSM = NULL
IBS_GEE = NULL
IBS_Gehan = NULL
IBS_SN = NULL

for(i in 1 : 10){
	Surv_obj = Surv(test_data[[i]]$Survival.Months, test_data[[i]]$Status=="Dead")

	pred_normal = NULL
	pred_NSNSM = NULL
	pred_Gehan = NULL
	pred_GEE = NULL
	pred_SN = NULL

	IBSrange = seq(from = 1, to = 90, length = 10)
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





########################################################