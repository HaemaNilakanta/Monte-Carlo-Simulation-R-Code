## Simple importance sampler example R code ##

	# set seed for reproducibility 
	set.seed(511)

	# Monte Carlo sample size
	m = 1e5

	# Classical Monte Carlo Estimate
	mean((1/3)*(rnorm(m, mean=0, sd=4) + rexp(m, 1) + rgamma(m, 1, 1/5)))

	# true mixture density
	f.true.density <- function(x){
	    return((1/3)*(dnorm(x, mean=0, sd=4) + dexp(x, 1) + dgamma(x, 1, 1/5)))
	}

	# draws from student's t distribution 
	x.vals = rt(m, df=1)

	# SIS estimate
	(muhat = mean(x.vals*f.true.density(x.vals)/dt(x.vals, df=1)))

	# sigma hat
	(sigmahat = sd(x.vals*f.true.density(x.vals)/dt(x.vals, df=1)))

	# 90% confidence interval 
	muhat + c(-1,1)*qnorm(0.95)*sigmahat/sqrt(m)
