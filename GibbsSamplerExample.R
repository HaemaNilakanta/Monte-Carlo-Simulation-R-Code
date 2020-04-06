## Gibbs sampler example R code ##

## Gibbs sampler example R code ##

	# set seed for reproducibility 
	set.seed(631)

	y.vals = c(5,1,5,14,3,19,1,1,4,22)
	hospital.times = c(2.5, 0.5, 1.3, 3, 1, 3.1, 0.1, 0.5, 1.5, 6.2)

	n = length(y.vals)

	# hyper parameters 
	a = c = 1
	d = 0.1

	# Monte Carlo sample size
	m = 1e5

	# where to store draws 
	lambda.mat = matrix(0, nrow=m, ncol=n)
	beta.vals = numeric(m)

	# starting values
	lambda.mat[1,] = rep(1, n)
	beta.vals[1] = 1

		# gibbs sampler 
	for(t in 2:m){
	    lambda.mat[t,] = rgamma(n, shape = (a + y.vals), 
	                               rate = (hospital.times + beta.vals[t-1]))
	    beta.vals[t] = rgamma(1, shape = (n*a + c), rate = (d + sum(lambda.mat[t-1,])))
	}

	# posterior means 
	c(colMeans(lambda.mat), mean(beta.vals))
	
	# plot of subset
	par(mfrow=c(2,2))
	hist(lambda.mat[,1], main="", freq=FALSE)
	hist(lambda.mat[,2], main="", freq=FALSE)
	hist(lambda.mat[,3], main="", freq=FALSE)
	hist(beta.vals, main="", freq=FALSE)
	par(mfrow=c(1,1))
