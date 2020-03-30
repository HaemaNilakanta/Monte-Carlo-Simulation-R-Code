## Linchpin Variable example R code ##

	# set seed for reproducibility 
	set.seed(331)

	# Generate data 
	n=100; p=3
	X = matrix(rnorm(n*p), nrow=n, ncol=p)
	beta.true = c(3.5, 1, -0.5)
	sigma.true = 1
	y.vals = X%*%beta.true + rnorm(n, mean=0, sd=sigma.true)

	# hyperparameters 
	a = b = 1

	# Posterior distribution terms
	A = t(X)%*%X + diag(p)
	A.inv = solve(A)
	A.svd = svd(A.inv) #singular value decomposition of A


	# Monte Carlo sample size
	m = 1e5

	# Draw from sigma^2|y
	sigma.rate.num = t(y.vals) %*% (diag(n) - X%*%A.inv%*%t(X)) %*% y.vals
	sigma.vals = 1 / rgamma(m, shape = (n/2 + a), rate = (sigma.rate.num/2 + b)) 

	# Draw from beta|sigma^2, y
	beta.mean = A.inv%*%t(X)%*%y.vals
	beta.varterm1 = A.svd$u %*% diag(sqrt(A.svd$d)) %*% t(A.svd$v)

	# store beta draws 
	beta.vals = matrix(0, ncol=p, nrow=m)
	for(t in 1:m){
	    beta.vals[t,] = beta.mean + beta.varterm %*% rnorm(p, sd=sqrt(sigma.vals[t]))
	}

	# estimated posterior means
	c(colMeans(beta.vals), mean(sigma.vals))

	# 90% credible intervals 
	beta90credible = apply(beta.vals, 2, function(x){quantile(x,c(0.1,0.9))})

	# append on credible interval for sigma. 
	cbind(beta90credible, quantile(sigma.vals, c(0.1,0.9))) 

# Plotting results 
	
	# histograms in a 2x2 grid with estimated density overlay
	par(mfrow=c(2,2))
	hist(beta.vals[,1], freq=FALSE, main="")
	    lines(density(beta.vals[,1]), lty=2)
	hist(beta.vals[,2], freq=FALSE, main="")
	    lines(density(beta.vals[,2]), lty=2)
	hist(beta.vals[,3], freq=FALSE, main="")
	    lines(density(beta.vals[,3]), lty=2)
	hist(sigma.vals, freq=FALSE, main="")
	    lines(density(sigma.vals), lty=2)
	par(mfrow=c(1,1))