## Gibbs sampler example R code ##

	# set seed for reproducibility 
	set.seed(631)

	# observed data
	y.vals = c(5,1,5,14,3,19,1,1,4,22)
	failuretime = c(94.320, 15.720, 62.880, 125.760,5.240,31.440, 1.048,1.048,2.096,10.480)
	n = length(y.vals)

	# hyper parameters 
	a = c = 1
	d = 0.1

	# Monte Carlo sample size
	m = 1e4

	# where to store draws 
	lambda.mat = matrix(0, nrow=m, ncol=n)
	beta.vals = numeric(m)

	# starting values
	lambda.mat[1,] = rep(1, n)
	beta.vals[1] = 1

	# gibbs sampler 
	for(t in 2:m){
	    lambda.mat[t,] = rgamma(n, shape = (a + y.vals), 
	                               rate = (failuretime + beta.vals[t-1]))
	    beta.vals[t] = rgamma(1, shape = (n*a + c), rate = (d + sum(lambda.mat[t-1,])))
	}

	# posterior means 
	c(colMeans(lambda.mat), mean(beta.vals))