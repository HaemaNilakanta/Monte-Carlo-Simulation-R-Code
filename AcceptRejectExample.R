## Accept-Reject example R code ##

# Generate draws from a truncated Normal(20,4) disribution where the values are greater than 28

# X ~ N(20, 4) then X = 20 + 4Z where Z ~ N(0,1)
	mean_hours = 20
	sd_hours = 4

	#Plot of N(20,4) with upper tail >28 shaded 
	x=seq(0, 40, length=400)
	y=dnorm(x, mean=mean_hours, sd = sd_hours)
	plot(x,y,type="l", lwd=2, col="black", xlab="Hours Studied Per Week", ylab="", yaxt='n')
	x.subset = seq(28.01, 39.99, length=200)
	y.subset = dnorm(x.subset, mean=mean_hours, sd = sd_hours)
	polygon(c(28,x.subset,40),c(0,y.subset,0),col="gray")


# Sample from a standard Normal distribution, only the upper tail when X>2

# set seed for reproducibility 
	set.seed(321)

	trunc.point = 2
	M = exp(-0.5*trunc.point^2)/(trunc.point*sqrt(2*pi)*pnorm(trunc.point, lower.tail=FALSE))

	#mean acceptance rate 
	print(1/M)

	#Monte Carlo sample size
	m = 1e4

	#truncated normal density
	f.truncatednormal <- function(x,c, mu, sigma){
	     return((dnorm(x, mean=mu, sd=sigma)/(1-pnorm(c, mean=mu, sd=sigma)))*ifelse(x>=c, 1, 0))
	}

	# generate iid realizations from Exp(1)
	xstar.vals1 = -log(1-runif(m))

	xstar.vals = trunc.point + xstar.vals1/trunc.point  #shift and scale by truncation point

	# independently generate iid realizations from U(0,1)
	u.vals = runif(m)

	# calculate the value of the densities 
	fstar.vals = trunc.point*exp(-trunc.point*(xstar.vals-trunc.point))
	ftruncated.vals = f.truncatednormal(x=xstar.vals, c=trunc.point, mu=0, sigma=1)

	acceptindices = u.vals < (ftruncated.vals/(M*fstar.vals))

	# proportion of acceptances
	mean(acceptindices)

	# keep only the accepted draws
	truncatedvals = xstar.vals[acceptindices]

	# put back in the original scale 
    truncatedvals_originalscale = mean_hours + truncatedvals*sd_hours

# Plotting results

    # Plotting on the standard normal scale 
   	# x-axis range for plotting and true density values
     xrange = seq(1.9, 5, by=0.01)
     f.true = f.truncatednormal(xrange,c=2, mu=0, sigma=1)

	#histogram of the output 
	hist(truncatedvals, freq=FALSE, breaks=30, 
	     xlim=c(1.9, 5), ylim=c(0,max(f.true)), 
	     main="")
	   lines(f.true~xrange, col="blue", lty=1, lw=3) #true density 

	# Plotting on the original scale 
	 xrange_originalscale = xrange*sd_hours + mean_hours
	 f.true_originalscale = f.truncatednormal(xrange_originalscale, c=28, mu=mean_hours, sigma=sd_hours)

	 	#histogram of the output 
	 pdf("histAcceptReject_originalscale.pdf")
	hist(truncatedvals_originalscale, freq=FALSE, breaks=30, 
	     xlim=c(min(xrange_originalscale), max(xrange_originalscale)), ylim=c(0,max(f.true_originalscale)), 
	     main="", xlab="Hours Studied Per Week")
	   lines(f.true_originalscale~xrange_originalscale, col="blue", lty=1, lw=3) #true density 
	   dev.off()

# Monte Carlo estimate and standard error #

	# Standard normal results 
		# MC estimate 
		# mean of truncated normal
		mean(truncatedvals)

		# MCSE of estimated mean
		mcse_tn = sd(truncatedvals)/sqrt(m)

		# 90% confidence interval, round to third decimal point
		round(mean(truncatedvals)+c(-1,1)*qnorm(0.95)*mcse_tn,3)

	# Original scale results
		# MC estimate 
		# mean of truncated normal
		mean(truncatedvals_originalscale)

		# MCSE of estimated mean
		mcse_tn_originalscale = sd(truncatedvals_originalscale)/sqrt(m)

		# 90% confidence interval, round to third decimal point
		round(mean(truncatedvals_originalscale)+c(-1,1)*qnorm(0.95)*mcse_tn_originalscale,3)
