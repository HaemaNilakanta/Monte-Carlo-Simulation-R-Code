## Weighted importance sampler example R code ##

	# set seed for reproducibility 
	set.seed(521)

	n = 50; p = 2

	# Monte Carlo sample size
	m = 1e4

	# simulate data 
	X = matrix(rnorm(n*p), nrow=n, ncol=p)
	beta.true = c(2, -1)

	# logit function
	prob.vec.true = exp(X%*%beta.true)/(1+exp(X%*%beta.true))

	# response variable 
	y.vec = rbinom(n, size=1, prob=prob.vec.true)

	glm.out = glm(y.vec ~ -1 + X, family="binomial")
	glm.out$coefficients

	# create a function for easy logit calculations 
	logit.func <- function(y, x, beta){
		return(exp(y*beta%*%x)/(1+exp(beta%*%x)))
	}

	# unnormalized posterior 
	log.hpost <- function(out, beta){
		unnormalized.posterior.log = (-sum(beta^2)/2) + sum(log(out))
		return(unnormalized.posterior.log)
	}

	#draws from the importance distribution 
	beta.draws = matrix(rnorm(m*p, mean=glm.out$coefficients, sd=1), 
	                              nrow = m, ncol=p, byrow=TRUE)
	#calculate the log-density
	log.ftilde.beta = apply(apply(beta.draws, 1, 
	                    function(x){log(dnorm(x, mean=glm.out$coefficients, sd=1))})
	                    ,2, sum)

	# store values
	out = matrix(NA, nrow=m, ncol=n)
	log.hpost.vec = numeric(m)

	for(t in 1:m){
	    for(i in 1:n){
	       out[t,i] = logit.func(y=y.vec[i], x=X[i,], beta=beta.draws[t,])
	    }
	   log.hpost.vec[t] = log.hpost(out = out[t,], beta = beta.draws[t,])
	}

	# weights 
	imp_weights = exp(log.hpost.vec - log.ftilde.beta)
	sum_of_imp_weights = sum(imp_weights)
	normalized_weights = imp_weights/sum_of_imp_weights

	# ratios for beta_1
	beta1_wis = beta.draws[,1]*imp_weights/sum_of_imp_weights
	beta1_squared_wis = (beta.draws[,1]^2)*imp_weights/sum_of_imp_weights

	# ratios for beta_2
	beta2_wis = beta.draws[,2]*imp_weights/sum_of_imp_weights
	beta2_squared_wis = (beta.draws[,2]^2)*imp_weights/sum_of_imp_weights

	#compute first and second moment of beta_2 with Monte Carlo standard errors
	sum(beta1_wis)
	sqrt(sum((normalized_weights^2)*(beta.draws[,1]-sum(beta1_wis))^2))

	sum(beta1_squared_wis)
	sqrt(sum((normalized_weights^2)*(beta.draws[,1]^2-sum(beta1_squared_wis))^2))

	#compute first and second moment of beta_1 with Monte Carlo standard errors
	sum(beta2_wis)
	sqrt(sum((normalized_weights^2)*(beta.draws[,2]-sum(beta2_wis))^2))

	#compute first and second moment of beta_2 with Monte Carlo standard errors
	sum(beta2_squared_wis)
	sqrt(sum((normalized_weights^2)*(beta.draws[,2]^2-sum(beta2_squared_wis))^2))