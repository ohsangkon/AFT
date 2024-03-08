###########################################################
## (Function) AFT with semiparametric skew normal scale mixture error  ##
###########################################################
###  Main function
#############################################
library(survival)
library(tidyverse)
library(sn)
library(nnls)

# Function

## 1. Log-likelihood of censored NSMSN
log_likelihood_AFTNSMSN = function(Eps, delta, Q_support, small_c)
{
  likelihood = density_AFTNSMSN(Eps = Eps, 
                                delta = delta,
                                Q_support = Q_support,
                                small_c = small_c)
  
  log_likelihood = likelihood %>% log() %>% sum()
  
  return(log_likelihood)
}

## 2. Density of NSMSN
density_AFTNSMSN = function(Eps, delta, Q_support, small_c)
{
  Eps_complete = Eps[which(delta == 1)]
  Eps_censored = Eps[which(delta != 1)]
  
  likelihood = rep(0, times = length(Eps))

  # xi's and lambda's are the same of all components, respectively.
  xi = Q_support[[1]]$xi
  lambda = Q_support[[1]]$lambda
  
  for (k in 1:length(Q_support))
  {
    prob_k = Q_support[[k]]$prob
    sigma_k = Q_support[[k]]$sigma + small_c
    
    if (length(Eps_complete) != 0)
    {
      likelihood_complete = dsn(x = Eps_complete, 
                                xi = xi,
                                omega = sigma_k,
                                alpha = lambda)*prob_k
      
      likelihood[which(delta == 1)] = 
        likelihood[which(delta == 1)] + likelihood_complete
    }
    if (length(Eps_censored) != 0)
    {
      likelihood_censored = (1 - psn(x = Eps_censored,
                                     xi = xi,
                                     omega = sigma_k,
                                     alpha = lambda))*prob_k
      
      likelihood[which(delta != 1)] =
        likelihood[which(delta != 1)] + likelihood_censored 
    }
  }
  
  return(likelihood)
}

## 3. Directional Derivative for censored SMSN
DD_AFTNSMSN = function(sigma, Eps, delta, Q_support, small_c)
{
  DD_denominator = 
    density_AFTNSMSN(Eps = Eps,
                     delta = delta,
                     Q_support = Q_support,
                     small_c = small_c)
  
  Q_support_numerator = list()
  Q_support_numerator[[1]] = 
    list(prob = 1, 
         xi = Q_support[[1]]$xi,
         sigma = sigma + small_c,
         lambda = Q_support[[1]]$lambda)
  
  DD_numerator = 
    density_AFTNSMSN(Eps = Eps,
                     delta = delta,
                     Q_support = Q_support_numerator,
                     small_c = small_c)
  n = length(DD_numerator)
  DD_value = sum(DD_numerator/DD_denominator) - n
  
  return(DD_value)
}

## 4. Empirical Derivative with central difference method
central_difference = function(func, value, h = 1e-4, ...)
{
  if ((length(value) != 1)) 
  {
    print("value should be single numeric")
    return(NULL)
  }
  numerator = func(value + h, ...) - func(value - h, ...)
  denominator = 2*h
  
  result = numerator/denominator
  
  return(result)
}

## 5. Find the grids that have local maximum.
where_potential_maximum_grid = function(logical_vector)
{
  target_vec = logical_vector
  len_vec = length(target_vec)
  
  A_vec = target_vec[-len_vec]
  B_vec = target_vec[-1]
  
  C_vec = (A_vec == T) & (B_vec == F)
  
  start_ind_2 = (1:(len_vec-1))[C_vec]
  end_ind_2 = start_ind_2 + 1
  
  if (target_vec[1] == F) {
    start_ind_1 = 1
    end_ind_1 = 2
  } else {
    start_ind_1 = NULL
    end_ind_1 = NULL
  }
  
  if (target_vec[len_vec] == T) {
    start_ind_3 = len_vec - 1
    end_ind_3 = len_vec
  } else {
    start_ind_3 = NULL
    end_ind_3 = NULL
  }
  
  start_ind = c(start_ind_1, start_ind_2, start_ind_3)
  end_ind = c(end_ind_1, end_ind_2, end_ind_3)
  
  mat_ind = cbind(start_ind, end_ind)
  colnames(mat_ind) = c("start", "end")
  
  return(mat_ind)
}

## 6. Make S Matrix
S_AFTNSMSN = function(Eps, delta, Q_support, small_c, sigma_new)
{
  S_denominator = 
    density_AFTNSMSN(Eps = Eps, 
                     delta = delta, 
                     Q_support = Q_support,
                     small_c = small_c)
  
  S_denominator_mat = 
    matrix(data = S_denominator,
           nrow = length(Eps),
           ncol = length(Q_support) + length(sigma_new))
  
  S_numerator_mat = 
    matrix(data = 0,
           nrow = length(Eps),
           ncol = length(Q_support) + length(sigma_new))
  
  Eps_complete = Eps[which(delta == 1)]
  Eps_censored = Eps[which(delta != 1)]
  
  max_k = length(Q_support) + length(sigma_new)
  xi_k = Q_support[[1]]$xi
  lambda_k = Q_support[[1]]$lambda
  
  for (k in 1:max_k)
  {
    S_k_numerator = rep(0, times = length(Eps))
    
    if (k <= length(Q_support))
    {
      sigma_k = Q_support[[k]]$sigma + small_c
    } else {
      sigma_k = sigma_new[(k - length(Q_support))] + small_c
    }
    
    if (length(Eps_complete) != 0)
    {
      S_k_complete = dsn(x = Eps_complete, 
                         xi = xi_k,
                         omega = sigma_k,
                         alpha = lambda_k)
      
      S_k_numerator[which(delta == 1)] = 
        S_k_numerator[which(delta == 1)] + S_k_complete 
    }
    
    if (length(Eps_censored) != 0)
    {
      S_k_censored = (1 - psn(x = Eps_censored,
                              xi = xi_k,
                              omega = sigma_k,
                              alpha = lambda_k))
      
      S_k_numerator[which(delta != 1)] =
        S_k_numerator[which(delta != 1)] + S_k_censored  
    }
    
    S_numerator_mat[,k] = S_k_numerator
  }
  
  S_mat = S_numerator_mat/S_denominator_mat
  
  return(S_mat)
}

## 7. Log_likelihood with Regression
## This is for an objective function for the function 'optim'
# beta_lambda is the vector that can be understood as below,
# beta_lambda = c(beta, Q_support[[1]]$lambda).

log_lik_reg_AFTNSMSN = function(X, beta_lambda, y, delta, Q_support, small_c)
{
  beta_temp = beta_lambda[1:ncol(X)]
  lambda_temp = beta_lambda[ncol(X) + 1]

  max_k = length(Q_support)
  mean_Q = 0
  for (k in 1:max_k)
  {
    prob_k = Q_support[[k]]$prob
    sigma_k = Q_support[[k]]$sigma + small_c
    
  }

 
  Q_support_temp = Q_support
  
  for (k in 1:length(Q_support))
 {
    Q_support_temp[[k]]$lambda = lambda_temp

  }
 
  Eps_temp = y - X %*% beta_temp
  
  log_lik = 
    log_likelihood_AFTNSMSN(Eps = Eps_temp,
                            delta = delta,
                            Q_support = Q_support_temp,
                            small_c = small_c)
  
  return(log_lik)
}


mean_Q_support = function(Q_support, small_c)
{
  max_k = length(Q_support)
  mean_Q = 0
  for (k in 1:max_k)
  {
    prob_k = Q_support[[k]]$prob
    sigma_k = Q_support[[k]]$sigma + small_c
    lambda = Q_support[[k]]$lambda
    
    mean_Q = mean_Q + prob_k*((sqrt(2/pi)*sigma_k*lambda/sqrt(1+lambda^2)))
  }
  
  return(mean_Q)
}


##########################################################################################
### Main Function  #############
##########################################################################################

AFT_NSNSM = function(x, y, delta = NULL, beta_hat = NULL, sigma_hat = NULL, lambda_hat = NULL, small_c = NULL, grid_max = NULL, largest_gradient_threshold = 0.1, tol = 1e-4, simulation_index = NULL){
  ### AFT using Nonparametric skew normal scale mixture error###
	library(tidyverse)
	library(sn)
	library(e1071)
	library(survival)

	if(is.null(delta) == TRUE){
		delta = matrix(1, nrow = nrow(x), ncol = 1)
	}

   #### Initial Parameters ####
	# beta & sigma & lambda   
	aft_sn = AFT_SN(x = x, y = y, delta = delta)  

	if(is.null(beta_hat) == TRUE){
		beta_hat = aft_sn[[1]] #aft_sn[[1]][-1]
	}

	if(is.null(sigma_hat) == TRUE){
		sigma_hat = aft_sn[[2]]
	}

	if(is.null(lambda_hat) == TRUE){
		lambda_hat = aft_sn[[3]]/2
	}

	if(is.null(small_c) == TRUE){
		small_c = aft_sn[[2]]/7
	}

	Q_support = list()
	Q_support[[1]] = list(prob = 1, xi = 0, sigma = sigma_hat, lambda = lambda_hat)

	# xi  
	xi_hat = -mean_Q_support(Q_support, 0)

	beta_hat[1] = beta_hat[1] + xi_hat

	#log likelihood
  	X = cbind(1,x)
	Eps = y - (X %*% beta_hat)

  log_likelihood = NULL
  log_likelihood[1] = 
    log_likelihood_AFTNSMSN(Eps = Eps, 
                            delta = delta,
                            Q_support = Q_support,
                            small_c = small_c)
  
  log_lik_diff = 1000
  i_AFTNSMSN = 0
  log_lik_diff_cnm = 1000
	diff_BFGS = 1000

	### Numerical method ###
    grid_max =max(100*IQR(y),2*(diff(range(y)))^2)
	 grid=seq(0,grid_max,(grid_max-0)/100)
 	 grid=c(grid,seq(grid[1],grid[2],(grid[2]-grid[1])/20))
	 sigma = NULL
	 for(i in 1:length(Q_support)){
		sigma = append(sigma, Q_support[[i]]$sigma)
	 }
 	 grid=sort(unique(c(0,grid,sigma)) )

  

while ((log_lik_diff_cnm > tol) | (diff_BFGS > tol))
{
  i_AFTNSMSN = i_AFTNSMSN + 1
  cat("##########  ", i_AFTNSMSN, "th AFTSMSN starts  ##########",
      "\n", sep = "")

	if(i_AFTNSMSN > 50 & log_lik_diff < tol){
		break
	}

  # CNM
  i_cnm = 0
  largest_gradient_value = 10
  Q_support_new = Q_support

  while ((i_cnm < 30) & (largest_gradient_value > largest_gradient_threshold))
  {
    i_cnm = i_cnm + 1
    cat(i_cnm, "th CNM \n", sep = "")

	 grid_for_DD = grid
    local_max_set = NULL
    DD_value_set = NULL
    
    # Originally, CNM use Newton + Bisection
    # For here, we just use the function "optimize"
    
    derivatives_of_gradient_function =
      sapply(X = grid_for_DD, 
             FUN = function(value){
               central_difference(func = DD_AFTNSMSN, value = value,
                                  Eps = Eps, Q_support = Q_support_new, small_c = small_c,
                                  delta = delta)
             })

	if(derivatives_of_gradient_function[length(derivatives_of_gradient_function)] > 0){

		 grid_max = grid_max + 1
		 grid=seq(0,grid_max,(grid_max-0)/100)
		 grid=c(grid,seq(grid[1],grid[2],(grid[2]-grid[1])/20))
		 sigma = NULL
		 for(i in 1:length(Q_support)){
			sigma = append(sigma, Q_support[[i]]$sigma)
		 }
	 	 grid=sort(unique(c(0,grid,sigma)) )

		 cat("Maximum of scale range:", " ", grid_max, "\n")
	}

    derivatives_positive = derivatives_of_gradient_function > 0
    
    grids_with_potential_mat = where_potential_maximum_grid(logical_vector = derivatives_positive)
    
    gradient_results = 
      apply(X = grids_with_potential_mat, MARGIN = 1,
            FUN = function(index){
              optimize(f = DD_AFTNSMSN, maximum = TRUE,
                       Eps = Eps, Q_support = Q_support_new, small_c = small_c,
                       delta = delta,
                       interval = c(grid_for_DD[index[1]], grid_for_DD[index[2]]))
            })
    
    gradient_results_mat = matrix(unlist(gradient_results), ncol = 2, byrow = T)
    colnames(gradient_results_mat) = c("maximum", "objective")
    
    index_maximum = which.max(gradient_results_mat[,"objective"])
    
    largest_gradient_value = gradient_results_mat[index_maximum, "objective"]
    if(largest_gradient_value < largest_gradient_threshold)
    {
      print("No More CNM iterations due to small gradient value")

		if(i_cnm == 1){
			    log_likelihood_new_cnm = log_likelihood_AFTNSMSN(Eps = Eps, 
                                                 delta = delta,
                                                 Q_support = Q_support_new,
                                                 small_c = small_c)
		}

      break
    }
      
    
    sigma_new = gradient_results_mat[, "maximum"]
    
    S_mat = 
      S_AFTNSMSN(Eps = Eps,
                 delta = delta,
                 Q_support = Q_support_new, small_c = small_c,
                 sigma_new = sigma_new)
    
    gamma_tol = 1e-10
    
    x_4nnls = rbind(S_mat*sqrt(gamma_tol), rep(1, ncol(S_mat)))
    
    y_4nnls = c(rep(2,nrow(S_mat))*sqrt(gamma_tol), 1)
    
    nnls_result = nnls(x_4nnls, y_4nnls)
    
    Q_support_new2 = Q_support_new
    
    nnls_result$x = (nnls_result$x)/sum(nnls_result$x)
    
    for (k in 1:(length(Q_support_new) + length(sigma_new)))
    {
      if (k <= length(Q_support_new)) 
      {
        Q_support_new2[[k]]$prob = nnls_result$x[k]
      } else {
        Q_support_new2[[k]] = list()
        Q_support_new2[[k]]$prob = nnls_result$x[k]
        Q_support_new2[[k]]$xi = Q_support_new[[1]]$xi
        Q_support_new2[[k]]$sigma = sigma_new[(k - length(Q_support_new))]
        Q_support_new2[[k]]$lambda = Q_support_new[[1]]$lambda
      }
    }
    
    idx_zero = (nnls_result$x == 0)
    
    Q_support_new2[idx_zero] = NULL
    
    log_likelihood_old_cnm = log_likelihood_AFTNSMSN(Eps = Eps, 
                                                 delta = delta,
                                                 Q_support = Q_support_new,
                                                 small_c = small_c)

    log_likelihood_new_cnm = log_likelihood_AFTNSMSN(Eps = Eps, 
                                                 delta = delta,
                                                 Q_support = Q_support_new2,
                                                 small_c = small_c)
    
    log_lik_diff_cnm = log_likelihood_new_cnm - log_likelihood_old_cnm

    num_support_new = length(Q_support_new2)
    cat("Largest DD: ", largest_gradient_value, 
        " Log-lik old: ", log_likelihood_old_cnm,
        " Log-lik new: ", log_likelihood_new_cnm,
        " Num_support_new: ", num_support_new, "\n", sep = "")

    if ((log_lik_diff_cnm) < 0)
    {
   	print("Something Wrong") 
      print("Log-likelihood usually increases monotonically")

		log_likelihood_new_cnm = log_likelihood_old_cnm

		break
    }

    Q_support_new = Q_support_new2

    if ((log_lik_diff_cnm) < tol) #$1e-10)
    {
	   	print("Loglikelihood is not increased signicantly") 
			break   
    }
      

  }   # CNM is terminated

  log_lik_diff_cnm = log_likelihood_new_cnm - log_likelihood[length(log_likelihood)]
  
	if(i_cnm == 1){
		  cat("# By the CNM, the log-likelihood increased as much as: ", 0, "\n", sep = "")
	}else{
	  cat("# By the CNM, the log-likelihood increased as much as: ", log_lik_diff_cnm, "\n", sep = "")
  }

  beta_lambda_temp = c(beta_hat, Q_support_new[[1]]$lambda)

  log_likelihood_old_reg =
    log_lik_reg_AFTNSMSN(X = X, beta_lambda = beta_lambda_temp,
                         y = y, delta = delta, Q_support = Q_support_new,
                         small_c = small_c)


  # Update beta & lambda 
     
  optim_result = 
    optim(par = beta_lambda_temp,
          fn = log_lik_reg_AFTNSMSN,
          X = X, y = y, delta = delta, 
          Q_support = Q_support_new, small_c = small_c,
          method = "BFGS", 
          control = list("fnscale" = -1))
  
  beta_new = optim_result$par[1:ncol(X)]
  lambda_new = optim_result$par[ncol(X) + 1]

  log_likelihood_new_reg = optim_result$value

  for (k in 1:length(Q_support_new))
  {
    Q_support_new[[k]]$lambda = lambda_new
  }

  Eps = y - (X %*% beta_new)

   log_likelihood_new_reg = log_likelihood_AFTNSMSN(Eps = Eps, 
                                                delta = delta,
                                                 Q_support = Q_support_new,
                                                 small_c = small_c)
  
  cat("Beta & Lambda Updated ", 
      " Log-lik old: ", log_likelihood_old_reg,
      " Log-lik new: ", log_likelihood_new_reg, "\n", sep = "")
 
	diff_BFGS = log_likelihood_new_reg - log_likelihood_old_reg

  cat("# By the BFGS, the log-likelihood increased as much as: ", diff_BFGS, "\n", sep = "")

	if(diff_BFGS < 0){
		if(abs(diff_BFGS) < tol ){
		  cat("It is may be occured by floating point: stop!", 
 				"\n", sep = "")
			converge = TRUE

		}else{
		  cat("BFGS algorithm may give another local solution: stop!", 
 				"\n", sep = "")
			converge = FALSE
		}
		break
	}

	log_likelihood = append(log_likelihood, log_likelihood_new_reg)

	log_lik_diff = log_likelihood[length(log_likelihood)] - log_likelihood[length(log_likelihood)-1]
  
		beta_hat = beta_new
		Q_support = Q_support_new

	if(log_lik_diff_cnm < tol | diff_BFGS < tol){
		converge = TRUE
	}

  } # Exit algorithm 

  #Sigma + small_c
  Q_support_new = Q_support

  for (k in 1:length(Q_support_new))
  {
    Q_support_new[[k]]$sigma = Q_support_new[[k]]$sigma + small_c
  }

  Q_support = Q_support_new


	# beta_0 & xi (adjustment)
	beta0_hat = beta_hat[1] + mean_Q_support(Q_support, 0)  
 	beta_hat = c(beta0_hat, beta_hat[-1])

  xi_new = -mean_Q_support(Q_support, 0)

  for (k in 1:length(Q_support))
  {
    Q_support[[k]]$xi = xi_new
  }


  result = list()
  result[[1]] = beta_hat
  result[[2]] = Q_support
  result[[3]] = log_likelihood
  result[[4]] = converge

  names(result) = c("beta_hat", "Other_estimates", "log_likelihood", "Converge")
  return(result) 
}


#########################################################################################

