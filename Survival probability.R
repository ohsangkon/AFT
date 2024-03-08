###################################################################################
# Compute P(survival time>t) given covariate
##################################################################################
library(sn)

survival_prob_NSNSM = function(t_star, model, covariate)
{

  #covariate: covariate for one observation
  #est_beta: estimated beta
  #lambda: estimated skewness

	est_beta = model$beta_hat

	Q_support = model$Other_estimates

	lambda = model$Other_estimates[[1]]$lambda

	if(is.null(nrow(covariate)) == FALSE ){

	  xx= as.matrix(cbind(1,covariate))
   	Eps_censored <<- log(t_star) - xx %*% est_beta

	}else{
		
		xx= c(1,covariate)
   	Eps_censored <<- log(t_star) - sum(xx * est_beta)
	}

  prob=0

  for (k in 1:length(Q_support))
   {
    prob = prob + (1 - psn(x = Eps_censored,
                                     xi = Q_support[[1]]$xi,
                                     omega = Q_support[[k]]$sigma,
                                     alpha = lambda))*Q_support[[k]]$prob
      
   }

  prob = as.matrix(prob)
  return(prob)
}

survival_prob_SN = function(t_star, model, covariate)
{
  #covariate: covariate for one observation
  #est_beta: estimated beta
  #lambda: estimated skewness

	est_beta = model[[1]] #model$beta_hat

	est_sigma = model[[2]]

	lambda = model[[3]]

	xi = model[[4]]

	if(is.null(nrow(covariate)) == FALSE ){

	  xx= as.matrix(cbind(1,covariate))
   	Eps_censored <<- log(t_star) - xx %*% est_beta

	}else{
		
		xx= c(1,covariate)
   	Eps_censored <<- log(t_star) - sum(xx * est_beta)
	}

    prob = (1 - psn(x = Eps_censored,
                    xi = xi,
                    omega = est_sigma,
                    alpha = lambda))
      

  prob = as.matrix(prob)
  return(prob)
}

survival_prob_Normal = function(t_star, model, covariate){

  #est_beta: estimated beta
  #est_sigma: estimated sigma

	est_beta = model$coefficients
	est_sigma = model$scale


	if(is.null(nrow(covariate)) == FALSE ){

	  xx= as.matrix(cbind(1,covariate))
   	Eps_censored <<- log(t_star) - xx %*% est_beta

	}else{
		
		xx= c(1,covariate)
   	Eps_censored <<- log(t_star) - sum(xx * est_beta)
	}

	prob = (1 - pnorm(q = Eps_censored, mean = 0, sd = est_sigma))

   return(prob)
}

survival_prob_Gehan = function(t_star, model_beta, covariate, train_x, train_y, d ){

  #est_beta: estimated beta
  #est_sigma: estimated sigma

	est_beta = model_beta

	if(is.null(nrow(covariate)) == FALSE ){

	  xx= as.matrix(cbind(1,covariate))
   	Eps_censored <<- log(t_star) - xx %*% est_beta

	}else{
		
		xx= c(1,covariate)
   	Eps_censored <<- log(t_star) - sum(xx * est_beta)
	}

	#ecdf from train data

	train_eps= (log(train_y) - as.matrix(cbind(1,train_x)) %*% est_beta)
	train_dt = as.data.frame(cbind(train_eps,d))
	colnames(train_dt) = c("res", "status")
	km = survfit(Surv(res,status)~1, data = train_dt)

	#train_fn = ecdf(train_eps)
	#train_prob = 1-train_fn(train_eps)

	#survial probability of test data
	prob = NULL
	for(i in 1 : length(Eps_censored)){
		if(is.na(which(Eps_censored[i] < km$time)[1]) == TRUE ){
			prob = append(prob, 0)
		}else if(which(Eps_censored[i] < km$time)[1] > 1){
			prob = append(prob, km$surv[which(Eps_censored[i] < km$time)[1]-1])
		}else{
			prob = append(prob, 1)
		}
	}

	prob = matrix(prob, ncol = 1)

	#fn = train_fn(Eps_censored)	
	#prob = 1-train_fn(Eps_censored)

	return(prob)
}

survival_prob_GEE = function(t_star, model, covariate, train_x, train_y, d){

  #est_beta: estimated beta
  #est_sigma: estimated sigma

	est_beta = model$coef.res

	#test
	if(is.null(nrow(covariate)) == FALSE ){

	  xx= as.matrix(cbind(1,covariate))
   	Eps_censored <<- log(t_star) - xx %*% est_beta

	}else{
		
		xx= c(1,covariate)
   	Eps_censored <<- log(t_star) - sum(xx * est_beta)
	}

	#ecdf from train data

	train_eps= (log(train_y) - as.matrix(cbind(1,train_x)) %*% est_beta)
	train_dt = as.data.frame(cbind(train_eps,d))
	colnames(train_dt) = c("res", "status")
	km = survfit(Surv(res,status)~1, data = train_dt)

	#survial probability of test data
	prob = NULL
	for(i in 1 : length(Eps_censored)){
		if(is.na(which(Eps_censored[i] < km$time)[1]) == TRUE ){
			prob = append(prob, 0)
		}else if(which(Eps_censored[i] < km$time)[1] > 1){
			prob = append(prob, km$surv[which(Eps_censored[i] < km$time)[1]-1])
		}else{
			prob = append(prob, 1)
		}
	}

	prob = matrix(prob, ncol = 1)

   return(prob)
}




