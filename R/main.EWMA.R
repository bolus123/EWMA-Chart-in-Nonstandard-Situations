require(pracma)

####################################################################################################################################################
    #parts of getting the charting constant l
####################################################################################################################################################

c4.f <- function(nu) sqrt(2 / nu) * 1 / beta(nu / 2, 1 / 2) * sqrt(pi)             #c4.function


####################################################################################################################################################
    # Phase II Xbar chart. Please see Yao and Chakraborti (2018)
####################################################################################################################################################

cs.f <- function(lambda, ss){

	numer <- 1 - (1 - lambda) ^ ss
	denom <- lambda / (2 - lambda) * (1 - (1 - lambda) ^ (2 * ss))
	
	numer / denom

}

EWMA.CFAR.DP.intgrand <- function(U, V, L, lambda, mm, ss) {

	nu <- mm - 1

	disp <- L * sqrt(qchisq(V, nu)) / c4.f(nu) / sqrt(mm - 1)
	
	ct <- cs.f(lambda, ss) * qnorm(U) / sqrt(mm)
	
	1 - pnorm(ct + disp) + pnorm(ct - disp)

}

EWMA.CARL.DP.intgrand <- function(U, V, L, lambda, mm, ss) {

	#cat(m, '\n')

	1 / (EWMA.CFAR.DP.intgrand(U, V, L, lambda, mm, ss)) 
	
}

#EWMA.CARL.DP.intgrand(3, 0.2, 50, Inf, u, v) 

EWMA.CARL.DP.integral <- function(xmin, xmax, ymin, ymax, L, lambda, mm, ss, reltol = 1e-6){

	#cat(m, '\n')

	integral2(EWMA.CARL.DP.intgrand, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, 
		singular = TRUE, L = L, lambda = lambda, mm = mm, ss = ss, reltol = reltol)$Q

}

EWMA.get.cc.DP <- function(ARL0 = 370, interval = c(1, 5), xmin = 0, xmax = 1, 
	ymin = 0, ymax = 1, lambda = 0.2, mm = 100, ss = Inf, reltol = 1e-6){

		root.finding <- function(ARL0, xmin, xmax, ymin, ymax, L, lambda, mm, ss, reltol){
		
			ARL0 - EWMA.CARL.DP.integral(xmin, xmax, ymin, ymax, L, lambda, mm, ss, reltol)
		
		}
	
		uniroot(root.finding, interval = interval, ARL0 = ARL0, xmin = xmin, xmax = xmax,
			ymin = ymin, ymax = ymax, lambda = lambda, mm = mm, ss = ss, reltol = reltol)$root
	
}

#plot.vec <- rep(NA, 6)
#
#lambda <- 0.9
#
#plot.vec[1] <- EWMA.get.cc.DP(ARL0 = 370, interval = c(1, 2.57), xmin = 0, xmax = 1, 
#	ymin = 0, ymax = 1, lambda = lambda, mm = 10, ss = Inf)
#	
#plot.vec[2] <- EWMA.get.cc.DP(ARL0 = 370, interval = c(1, 3), xmin = 0, xmax = 1, 
#	ymin = 0, ymax = 1, lambda = lambda, mm = 20, ss = Inf)
#	
#plot.vec[3] <- EWMA.get.cc.DP(ARL0 = 370, interval = c(1, 4), xmin = 0, xmax = 1, 
#	ymin = 0, ymax = 1, lambda = lambda, mm = 300, ss = Inf)
#	
#plot.vec[4] <- EWMA.get.cc.DP(ARL0 = 370, interval = c(1, 4), xmin = 0, xmax = 1, 
#	ymin = 0, ymax = 1, lambda = lambda, mm = 100, ss = Inf)	
#	
#plot.vec[5] <- EWMA.get.cc.DP(ARL0 = 370, interval = c(1, 4), xmin = 0, xmax = 1, 
#	ymin = 0, ymax = 1, lambda = lambda, mm = 200, ss = Inf)	
#	
#plot.vec[6] <- EWMA.get.cc.DP(ARL0 = 370, interval = c(1, 4), xmin = 0, xmax = 1, 
#	ymin = 0, ymax = 1, lambda = lambda, mm = 300, ss = Inf)	


	
####################################################################################################################################################

EWMA.ARLin.performance <- function(L, lambda, mm, ss, sim = 10000){

	U <- runif(sim)
	V <- runif(sim)

	out <- EWMA.CARL.DP.intgrand(U, V, L, lambda, mm, ss)
	
	res <- list(mean = mean(out), sd = sd(out), 
		quantiles = quantile(out, probs = c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1))
	)
	
	res
	
}

#EWMA.ARLin.performance(L = 2.963359, lambda = 0.2, mm = 20, ss = Inf, sim = 1000000)

####################################################################################################################################################

EWMA.CARL.MC.Q.TRAD <- function(mu, sigma, L, lambda, ss, tt) {

	nn <- 2 * tt + 1

	Q <- matrix(NA, ncol = nn, nrow = nn)
	
	for (l in 1:nn){
	
		ll <- l - (tt + 1)
	
		for (k in 1:nn){
		
			kk <- k - (tt + 1)
		
			Q[l, k] <- pnorm(mu / sigma + ( 2 * lambda * tt + lambda + 2 * ll - 2 * (1 - lambda) * kk  + 1) / nn / lambda * L * 
				sqrt(lambda / (2 - lambda) * (1 - (1 - lambda) ^ (2 * ss)))) - 
				pnorm(mu / sigma + (2 * lambda * tt + lambda + 2 * ll - 2 * (1 - lambda) * kk - 1) / nn / lambda * L * 
				sqrt(lambda / (2 - lambda) * (1 - (1 - lambda) ^ (2 * ss))))
		
		}
	
	}
	
	Q

}



EWMA.CARL.MC.xi <- function(tt, mid = tt + 1) {

	nn <- 2 * tt + 1
	
	xi.vec <- rep(0, nn)
	
	xi.vec[mid] <- 1

	xi.vec

}

EWMA.CARL.MC.TRAD.prod <- function(mu, sigma, L, lambda, ss, tt){

	QQ <- EWMA.CARL.MC.Q.TRAD(mu, sigma, L, lambda, ss, tt) 
	xi.vec <- EWMA.CARL.MC.xi(tt)
	xi.vec <- matrix(xi.vec, nrow = 1, ncol = 2 * tt + 1)
	I.matrix <- diag(2 * tt + 1)
	
	one.vec <- matrix(1, ncol = 1, nrow = 2 * tt + 1)
	
	inv.matrix <- try(solve(I.matrix - QQ), silent = TRUE)
	
	if (class(inv.matrix) == 'try-error') {
	
		cat('try-error', '\n')
		inv.matrix <- MASS::ginv(I.matrix - QQ)
	
	}
	
	xi.vec %*% inv.matrix %*% one.vec

}

#EWMA.CARL.MC.TRAD.prod(mu = 0, sigma = 1, L = 3, lambda = 1, ss = Inf, tt = 1000)


EWMA.get.cc.MC.TRAD <- function(ARL0 = 370, interval = c(1, 5), mu = 0, sigma = 1, 
	lambda = 0.2, ss = Inf, tt = 3, reltol = 1e-6, tol = 1e-6){

		root.finding <- function(ARL0, mu, sigma, L, lambda, ss, tt, reltol){
		
			ARLin <- EWMA.CARL.MC.TRAD.prod(mu, sigma, L, lambda, ss, tt)
			
			cat('ARLin:', ARLin, '\n')
			
			ARL0 - ARLin
		
		}
	
		uniroot(root.finding, interval = interval, tol = tol, ARL0 = ARL0, mu = mu, sigma = sigma, 
			lambda = lambda, ss = ss, tt = tt, reltol = reltol)$root
	
}

#EWMA.get.cc.MC.TRAD(ARL0 = 370, interval = c(1, 5), mu = 0, sigma = 1, 
#	lambda = 0.2, ss = Inf, tt = 3, reltol = 1e-6, tol = 1e-6)


EWMA.CARL.MC.Q <- function(U, V, L, lambda, mm, ss, tt) {

	nn <- 2 * tt + 1

	Q <- matrix(NA, ncol = nn, nrow = nn)
	
	for (l in 1:nn){
	
		ll <- l - (tt + 1)
	
		for (k in 1:nn){
		
			kk <- k - (tt + 1)
		
			Q[l, k] <- pnorm(qnorm(U) / sqrt(mm) + (2 * kk - (1 - lambda) * 2 * ll + 1) / nn * 
				(L * sqrt(1 / lambda / (2 - lambda) * (1 - (1 - lambda)^(2 * ss))) * sqrt(qchisq(V, mm - 1)) / c4.f(mm - 1) / sqrt(mm - 1) )) - 
				pnorm(qnorm(U) / sqrt(mm) + (2 * kk - (1 - lambda) * 2 * ll - 1) / nn * 
				(L * sqrt(1 / lambda / (2 - lambda) * (1 - (1 - lambda)^(2 * ss))) * sqrt(qchisq(V, mm - 1)) / c4.f(mm - 1) / sqrt(mm - 1) )) 
		
		}
	
	}
	
	Q
	
}

EWMA.CARL.MC.xi <- function(tt, mid = tt + 1) {

	nn <- 2 * tt + 1
	
	xi.vec <- rep(0, nn)
	
	xi.vec[mid] <- 1

	xi.vec

}

EWMA.CARL.MC.integrand <- function(U, V, L, lambda, mm, ss, tt){

	QQ <- EWMA.CARL.MC.Q(U, V, L, lambda, mm, ss, tt)
	xi.vec <- EWMA.CARL.MC.xi(tt)
	xi.vec <- matrix(xi.vec, nrow = 1, ncol = 2 * tt + 1)
	I.matrix <- diag(2 * tt + 1)
	
	one.vec <- matrix(1, ncol = 1, nrow = 2 * tt + 1)
	
	inv.matrix <- try(solve(I.matrix - QQ), silent = TRUE)
	
	if (class(inv.matrix) == 'try-error') {
	
		cat('try-error', '\n')
		inv.matrix <- MASS::ginv(I.matrix - QQ)
	
	}
	
	xi.vec %*% inv.matrix %*% one.vec

}

#debug(EWMA.CARL.MC)

#EWMA.CARL.MC(0.5, 0.5, 3, 0.2, 100, Inf, 5)

EWMA.CARL.MC.integral <- function(xmin, xmax, ymin, ymax, L, lambda, mm, ss, tt = 3, reltol = 1e-6){

	#cat(m, '\n')

	integral2(EWMA.CARL.MC.integrand, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, 
		singular = TRUE, vectorized = FALSE, L = L, lambda = lambda, mm = mm, ss = ss, tt = tt, reltol = reltol)$Q

}



#debug(EWMA.CARL.MC.Q)

EWMA.get.cc.MC <- function(ARL0 = 370, interval = c(1, 5), xmin = 0, xmax = 1, 
	ymin = 0, ymax = 1, lambda = 0.2, mm = 100, ss = Inf, tt = 3, reltol = 1e-6, tol = 1e-6){

		root.finding <- function(ARL0, xmin, xmax, ymin, ymax, L, lambda, mm, ss, tt, reltol){
		
			ARLin <- EWMA.CARL.MC.integral(xmin, xmax, ymin, ymax, L, lambda, mm, ss, tt, reltol)
			
			cat('L:', L, ', ARLin:', ARLin, '\n')
		
			ARL0 - ARLin
		
		}
	
		uniroot(root.finding, interval = interval, tol = tol, ARL0 = ARL0, xmin = xmin, xmax = xmax,
			ymin = ymin, ymax = ymax, lambda = lambda, mm = mm, ss = ss, tt = tt, reltol = reltol)$root
	
}

#setMKLthreads(getMKLthreads() - 1)

EWMA.get.cc.MC(ARL0 = 370, interval = c(2.3, 3.5), xmin = 0, xmax = 1, 
		ymin = 0, ymax = 1, lambda = 0.8, mm = 100, ss = Inf, tt = 50, reltol = 1e-6)

####################################################################################################################################################

#EWMA.sim <- function(L, lambda, mm, z.sim = 1000, sim = 1000){
#	
#	RL.vec <- rep(NA, sim)
#	
#	for (ii in 1:sim){
#
#		PH1.X.bar <- rnorm(mm)
#	
#		mu.hat <- mean(PH1.X.bar)
#		sigma.hat <- sd(PH1.X.bar) / c4.f(mm - 1)
#	
#		RL <- 0
#		RL.done <- 0
#		
#		Z.vec <- rep(NA, z.sim)
#		
#		while (RL.done == 0){
#		
#			X.bar.vec <- rnorm(z.sim, mu.hat, sigma.hat)
#		
#			Z.0 <- ifelse(RL == 0, mu.hat, Z.vec[z.sim])
#		
#			for (jj in 1:z.sim) {
#
#				Z.vec[jj] <- lambda * X.bar.vec[jj] + (1 - lambda) * Z.0
#				
#				Z.0 <- Z.vec[jj]
#			
#			}
#			
#			LCL <- mu.hat - L * sigma.hat * sqrt(lambda / (2 - lambda))
#			UCL <- mu.hat + L * sigma.hat * sqrt(lambda / (2 - lambda))
#			
#			IC <- LCL <= Z.vec & Z.vec <= UCL
#			
#			if (length(which(IC)) == z.sim) {
#			
#				RL <- RL + z.sim
#			
#			} else {
#			
#				RL <- RL + which.min(IC)
#				RL.done <- 1
#			
#			}
#		
#		}
#		
#		RL.vec[ii] <- RL
#		
#	}
#	
#	res <- list(ARL = mean(RL.vec), SDRL = sd(RL.vec), RL = RL.vec)
#	
#	return(res)
#	
#}	
#
##debug(EWMA.sim)
#
#a <- EWMA.sim(2.983686, 0.8, 300, z.sim = 1000, sim = 10000)
#
#l.vec <- c(
#	2.369803, 2.681122, 2.980097, 2.939529, 2.970110, 2.980097
#)
#
#lambda <- c(0.9)
#lambda.vec <- rep(lambda, 6)
#
#mm.vec <- c(10, 20, 50, 100, 200, 300)
#
#b.matrix <- matrix(NA, ncol = 6, nrow = 100)
#
#for (j in 1:100){
#
#	b.vec <- rep(NA, 6)
#
#	for (i in 1:6){
#		b.vec[i] <- EWMA.sim(l.vec[i], lambda.vec[i], mm.vec[i], z.sim = 1000, sim = 10000)$ARL
#	}
#	
#	b.matrix[j, ] <- b.vec
#
#
#}
#