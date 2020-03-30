## Accept-Reject example R code ##

# Sample from a standard Normal distribution, only the upper tail when X>2

# set seed for reproducibility 
	set.seed(321)

	trunc.point = 2
	M = exp(-0.5*trunc.point^2)/(trunc.point*sqrt(2*pi)*pnorm(trunc.point, lower.tail=FALSE))

	#mean acceptance rate 
	print(1/M)
	[1] 0.8427385

	#Monte Carlo sample size
	m = 1e4

	#truncated normal density
	f.truncatednormal <- function(x,c){
	     return((dnorm(x)/(1-pnorm(c)))*ifelse(x>=c, 1, 0))
	}

	# generate iid realizations from Exp(1)
	xstar.vals1 = -log(1-runif(m))

	xstar.vals = trunc.point + xstar.vals1/trunc.point  #shift and scale by truncation point

	# independently generate iid realizations from U(0,1)
	u.vals = runif(m)

	# calculate the value of the densities 
	fstar.vals = trunc.point*exp(-trunc.point*(xstar.vals-trunc.point))
	ftruncated.vals = f.truncatednormal(x=xstar.vals, c=trunc.point)

	acceptindices = u.vals < (ftruncated.vals/(M*fstar.vals))

	# proportion of acceptances
	mean(acceptindices)

	# keep only the accepted draws
	truncatedvals = xstar.vals[acceptindices]

# Plotting results

   # x-axis range for plotting and true density values
     xrange = seq(1.9, 5, by=0.01)
     f.true = f.truncatednormal(xrange,c=2)
     
	#histogram of the output 
	hist(truncatedvals, freq=FALSE, breaks=30, 
	     xlim=c(1.9, 5), ylim=c(0,max(f.true)), 
	     main="")
	   lines(f.true~xrange, col="blue", lty=1, lw=3) #true density 