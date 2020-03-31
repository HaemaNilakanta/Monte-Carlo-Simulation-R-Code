	## Metropolis-Hasting example R code ##

	# set seed for reproducibility 
	set.seed(621)

	# Use same data as Linchpin Method
	# create a log-likelihood function
	loglike.func <- function(y, X, beta, sigmasquared, a, b){
	    n = length(y)
	    p = ncol(X)

	    a1 = sum((y - X%*%beta)^2)/2
	    a2 = sum((beta)^2)/2
	    a3 = b

	    out = -0.5*(n + p + 2*a + 2)*log(sigmasquared) + (-a1-a2-a3)/sigmasquared
	    return(out)
	}

	# Monte Carlo sample size 
	m = 1e5 

	#hyper parameters 
	a = b = 1

	# create shells to hold coefficients after each iteration
	beta.mat = matrix(0, nrow=m, ncol=p)
	sigmasquared.vec = numeric(m)

	#run an OLS for starting values, -1 for no intercept
	lm.out = lm(y.vals ~ -1 + X)

	#keep track of acceptances/rejection 
	accept.vec = numeric(m)

	# set initial values 
	beta.mat[1,] = lm.out$coefficients
	sigmasquared.vec[1] = sigma(lm.out)
	accept.vec[1] = NA


	for(t in 2:m){
		#propose a move
		proposal.move = c(beta.mat[t-1,], sigmasquared.vec[t-1]) + rnorm((p+1), mean=0, sd=0.1)

		# already know that sigma is a variance term so can't be negative 
		# therefore if proposed move is negative, reject 
		if(proposal.move[p+1]<=0){
		    beta.mat[t,] = beta.mat[t-1,]
		    sigmasquared.vec[t] = sigmasquared.vec[t-1]
		    accept.vec[t] = 0
		} else{
		    #otherwise check the ratio 
		    logratio = loglike.func(y=y.vals, X=X, beta=proposal.move[1:p],
		                    sigmasquared=proposal.move[p+1], a=a, b=b) - 
			            loglike.func(y=y.vals, X=X, beta=beta.mat[t-1,], 
			                sigmasquared=sigmasquared.vec[t-1], a=a, b=b)

		    u.val = runif(1)

		    if(u.val < exp(logratio)){
		        beta.mat[t,] = proposal.move[1:p]
		        sigmasquared.vec[t] = proposal.move[p+1]
		        accept.vec[t] = 1
		    } else{
		        beta.mat[t,] = beta.mat[t-1,]
		        sigmasquared.vec[t] = sigmasquared.vec[t-1]
		        accept.vec[t] = 0
            }
        }

    }

    # mean acceptance rate
	mean(accept.vec, na.rm=TRUE)

	# histograms in a 2x2 grid
	par(mfrow=c(2,2))
	   hist(beta.mat[,1], main="", freq=FALSE)
	   hist(beta.mat[,2], main="", freq=FALSE)
	   hist(beta.mat[,3], main="", freq=FALSE)
	   hist(sigmasquared.vec, main="", freq=FALSE)
	par(mfrow=c(1,1))

#Compute MCSEs

	# load library
	library(mcmcse)

	#posterior means and Monte Carlo standard error
	(mcse_out = mcse.mat(beta.mat))

	# Compute confidnece intervals 
		# get t quantile, df=number of batches-1
	   batchsize = floor(sqrt(nrow(beta.mat)))
	   batches = floor(nrow(beta.mat)/batchsize)
	   tstar = qt(0.95, batches-1)

	# 90% CI for beta_1
	mcse_out[1,1] + c(-1,1)*tstar*mcse_out[1,2]

	# 90% CI for beta_2
	mcse_out[2,1] + c(-1,1)*tstar*mcse_out[2,2]

	# 90% CI for beta_3
	mcse_out[3,1] + c(-1,1)*tstar*mcse_out[3,2]