# AFT (Adaptive Accelerated Failure Time modeling with a Semiparametric Skewed Error Distribution)
The accelerated failure time (AFT) model is widely used to analyze relationships between variables in the presence of censored observations. However, this model relies on some assumptions such as the error distribution, which can lead to biased or inefficient estimates if these assumptions are violated. In order to overcome this challenge, we propose a novel approach that incorporates a semiparametric skewnormal scale mixture distribution for the error term in the AFT model. By allowing for more flexibility and robustness, this approach reduces the risk of misspecification
and improves the accuracy of parameter estimation. 

(Oh, S., Lee, H., Kang, S., and Seo, B. (2024). Adaptive accelerated failure time modeling with a semiparametric skewed error distribution. arXiv:2402.02128.)


# AFT_NSNSM.R
Basic codes for Adaptive AFT modeling with a semiparametric skew-normal scale mixture error distribution. There are 9 functions in this file. AFT_NSNSM function is the main function, and others are required to save if you want to use AFT_NSNSM function.

# AFT_SN.R
Codes for AFT modeling with a skew-normal error distribution. There are 4 functions in this file. AFT_SN function is the main function, and others are required to save if you want to use AFT_SN function. Furthermore, AFT_SN function is used for initial values for AFT_NSNSM function.

# Survival probability.R
Computation for P(survival time > t) given covariates.

survival_prob_NSNSM: Survival probability based on AFT with a semiparametric skew-normal scale mixture error

survival_prob_SN: Survival probability based on AFT with a skew-normal error

survival_prob_Normal: Survival probability based on AFT with a normal error

survival_prob_Gehan: Survival probability based on semiparametric AFT (Gehan) 

survival_prob_GEE: Survival probability based on semiparametric AFT (GEE) 

# Simulation
Codes for simulation studies:

simulation1

simulation2

simulation3

# Lung_cancer.R
Lung cancer dataset analysis using AFT models.

# breast_cancer.R
breast cancer dataset analysis using AFT models.

Breast_Cancer.csv should be loaded.
