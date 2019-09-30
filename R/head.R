require(pracma)

####################################################################################################################################################
    #parts of getting the charting constant l
####################################################################################################################################################

c4.f <- function(nu) sqrt(2 / nu) * 1 / beta(nu / 2, 1 / 2) * sqrt(pi)             #c4.function


####################################################################################################################################################
    # Phase II Xbar chart. Please see Yao and Chakraborti (2018)
####################################################################################################################################################

#cs.f <- function(lambda, ss){
#
#	numer <- 1 - (1 - lambda) ^ ss
#	denom <- lambda / (2 - lambda) * (1 - (1 - lambda) ^ (2 * ss))
#	
#	numer / denom
#
#}
#
#EWMA.CFAR.DP.intgrand <- function(U, V, L, lambda, mm, ss) {
#
#	nu <- mm - 1
#
#	disp <- L * sqrt(qchisq(V, nu)) / c4.f(nu) / sqrt(mm - 1)
#	
#	ct <- cs.f(lambda, ss) * qnorm(U) / sqrt(mm)
#	
#	1 - pnorm(ct + disp) + pnorm(ct - disp)
#
#}
#
#EWMA.CARL.DP.intgrand <- function(U, V, L, lambda, mm, ss) {
#
#	#cat(m, '\n')
#
#	1 / (EWMA.CFAR.DP.intgrand(U, V, L, lambda, mm, ss)) 
#	
#}

#EWMA.CARL.DP.intgrand(3, 0.2, 50, Inf, u, v) 

#EWMA.CARL.DP.integral <- function(xmin, xmax, ymin, ymax, L, lambda, mm, ss, reltol = 1e-6){
#
#	#cat(m, '\n')
#
#	integral2(EWMA.CARL.DP.intgrand, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, 
#		singular = TRUE, L = L, lambda = lambda, mm = mm, ss = ss, reltol = reltol)$Q
#
#}
#
#EWMA.get.cc.DP <- function(ARL0 = 370, interval = c(1, 5), xmin = 0, xmax = 1, 
#	ymin = 0, ymax = 1, lambda = 0.2, mm = 100, ss = Inf, reltol = 1e-6){
#
#		root.finding <- function(ARL0, xmin, xmax, ymin, ymax, L, lambda, mm, ss, reltol){
#		
#			ARL0 - EWMA.CARL.DP.integral(xmin, xmax, ymin, ymax, L, lambda, mm, ss, reltol)
#		
#		}
#	
#		uniroot(root.finding, interval = interval, ARL0 = ARL0, xmin = xmin, xmax = xmax,
#			ymin = ymin, ymax = ymax, lambda = lambda, mm = mm, ss = ss, reltol = reltol)$root
#	
#}

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

#EWMA.ARLin.performance <- function(L, lambda, mm, ss, sim = 10000){
#
#	U <- runif(sim)
#	V <- runif(sim)
#
#	out <- EWMA.CARL.DP.intgrand(U, V, L, lambda, mm, ss)
#	
#	res <- list(mean = mean(out), sd = sd(out), 
#		quantiles = quantile(out, probs = c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1))
#	)
#	
#	res
#	
#}

#EWMA.ARLin.performance(L = 2.963359, lambda = 0.2, mm = 20, ss = Inf, sim = 1000000)

####################################################################################################################################################

EWMA.CARL.MC.Q.TRAD <- function(mu, sigma, L, lambda, ss, tt) {

	nn <- 2 * tt + 1

	Q <- matrix(NA, ncol = nn, nrow = nn)
	
	zz <- seq(-tt, tt, 1)
	
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
	
		#cat('try-error', '\n')
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


EWMA.CARL.MC.Q <- function(U, V, L, lambda, mm, ss, tt, delta = 0) {

	nn <- 2 * tt + 1

	Q <- matrix(NA, ncol = nn, nrow = nn)
	
	kk <- seq(-tt, tt, 1)
	
	theta <- 1 - (1 - lambda) ^ (2 * ss)
	
	tau <- L * sqrt(qchisq(V, mm - 1) / (mm - 1)) / c4.f(mm - 1) * 
		sqrt(lambda / (2 - lambda)) * sqrt(theta) / (2 * tt + 1) 
	
	cent <- qnorm(U) / sqrt(mm)
	
	for (l in 1:nn){
	
		ll <- l - (tt + 1)
	
		#for (k in 1:nn){
		
		#	kk <- k - (tt + 1)
		
			#Q[l, ] <- pnorm(qnorm(U) / sqrt(mm) + (2 * zz - (1 - lambda) * 2 * ll + 1) / nn * 
			#	(L * sqrt(1 / lambda / (2 - lambda) * (1 - (1 - lambda)^(2 * ss))) * sqrt(qchisq(V, mm - 1)) / c4.f(mm - 1) / sqrt(mm - 1) ) - delta) - 
			#	pnorm(qnorm(U) / sqrt(mm) + (2 * zz - (1 - lambda) * 2 * ll - 1) / nn * 
			#	(L * sqrt(1 / lambda / (2 - lambda) * (1 - (1 - lambda)^(2 * ss))) * sqrt(qchisq(V, mm - 1)) / c4.f(mm - 1) / sqrt(mm - 1) ) - delta) 

			Q[, l] <- pnorm(cent + (2 * ll - (1 - lambda) * 2 * kk + 1) * tau / lambda - delta) - 
				pnorm(cent + (2 * ll - (1 - lambda) * 2 * kk - 1) * tau / lambda - delta)
		
		#}
	
	}
	
	Q
	
}

EWMA.CARL.MC.xi <- function(tt, mid = tt + 1) {

	nn <- 2 * tt + 1
	
	xi.vec <- rep(0, nn)
	
	xi.vec[mid] <- 1

	xi.vec

}

EWMA.CARL.MC.integrand <- function(U, V, L, lambda, mm, ss, tt, delta = 0){

	QQ <- EWMA.CARL.MC.Q(U, V, L, lambda, mm, ss, tt, delta = delta)
	xi.vec <- EWMA.CARL.MC.xi(tt)
	xi.vec <- matrix(xi.vec, nrow = 1, ncol = 2 * tt + 1)
	I.matrix <- diag(2 * tt + 1)
	
	one.vec <- matrix(1, ncol = 1, nrow = 2 * tt + 1)
	
	inv.matrix <- try(solve(I.matrix - QQ), silent = TRUE)
	
	if (class(inv.matrix) == 'try-error') {
	
		#cat('try-error', '\n')
		inv.matrix <- MASS::ginv(I.matrix - QQ)
	
	}
	
	xi.vec %*% inv.matrix %*% one.vec

}

#debug(EWMA.CARL.MC)

#EWMA.CARL.MC(0.5, 0.5, 3, 0.2, 100, Inf, 5)

EWMA.CARL.MC.integral <- function(xmin, xmax, ymin, ymax, L, lambda, mm, ss, tt = 3, delta = 0, reltol = 1e-6){

	#cat(m, '\n')

	integral2(EWMA.CARL.MC.integrand, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, 
		singular = TRUE, vectorized = FALSE, L = L, lambda = lambda, mm = mm, ss = ss, tt = tt, delta = delta, reltol = reltol)$Q

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


####################################################################################################################################################


EWMA.CARL.Conditional.MC.integrand <- function(U, V, L, lambda, eplison, ARL0, mm, ss, tt){

	QQ <- EWMA.CARL.MC.Q(U, V, L, lambda, mm, ss, tt)
	xi.vec <- EWMA.CARL.MC.xi(tt)
	xi.vec <- matrix(xi.vec, nrow = 1, ncol = 2 * tt + 1)
	I.matrix <- diag(2 * tt + 1)
	
	one.vec <- matrix(1, ncol = 1, nrow = 2 * tt + 1)
	
	inv.matrix <- try(solve(I.matrix - QQ), silent = TRUE)
	
	if (class(inv.matrix) == 'try-error') {
	
		#cat('try-error', '\n')
		inv.matrix <- MASS::ginv(I.matrix - QQ)
	
	}
	
	ifelse(xi.vec %*% inv.matrix %*% one.vec >= (1 - eplison) * ARL0, 1, 0)
	
	#ifelse(QQ >= (1 - eplison) * ARL0, 1, 0)

}


EWMA.CARL.Conditional.MC.integral <- function(xmin, xmax, ymin, ymax, L, lambda, eplison, ARL0, mm, ss, tt = 3, reltol = 1e-6){

	#cat(m, '\n')

	integral2(EWMA.CARL.Conditional.MC.integrand, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, 
		singular = FALSE, vectorized = FALSE, L = L, lambda = lambda, eplison = eplison, ARL0 = ARL0, 
		mm = mm, ss = ss, tt = tt, reltol = reltol)$Q

}


EWMA.get.cc.Conditional.MC <- function(p0 = 0.05, interval = c(1, 5), xmin = 0, xmax = 1, 
	ymin = 0, ymax = 1, lambda = 0.2, eplison = 0.1, ARL0 = 370, mm = 100, ss = Inf, tt = 3, reltol = 1e-6, tol = 1e-6){

		root.finding <- function(p0, xmin, xmax, ymin, ymax, L, lambda, eplison, ARL0, mm, ss, tt, reltol){
		
			p <- EWMA.CARL.Conditional.MC.integral(xmin, xmax, ymin, ymax, 
				L, lambda, eplison, ARL0, mm, ss, tt, reltol)
			
			cat('L:', L, ', p:', 1 - p, '\n')
		
			1 - p0 - p
		
		}
	
		uniroot(root.finding, interval = interval, tol = tol, p0 = p0, xmin = xmin, xmax = xmax,
			ymin = ymin, ymax = ymax, lambda = lambda, eplison = eplison, ARL0 = ARL0, 
			mm = mm, ss = ss, tt = tt, reltol = reltol)$root
	
}
#
#EWMA.get.cc.Conditional.MC(p0 = 0.1, interval = c(2, 4), xmin = 0, xmax = 1, 
#	ymin = 0, ymax = 1, lambda = 0.2, eplison = 0, ARL0 = 370, mm = 50, ss = Inf, tt = 50, reltol = 1e-6, tol = 1e-6)
#	
#EWMA.get.cc.Conditional.MC(p0 = 0.1, interval = c(2, 4), xmin = 0, xmax = 1, 
#	ymin = 0, ymax = 1, lambda = 0.2, eplison = 0, ARL0 = 370, mm = 50, ss = Inf, tt = 100, reltol = 1e-6, tol = 1e-6)
#	
#EWMA.get.cc.Conditional.MC(p0 = 0.1, interval = c(2, 4), xmin = 0, xmax = 1, 
#	ymin = 0, ymax = 1, lambda = 0.2, eplison = 0, ARL0 = 370, mm = 50, ss = Inf, tt = 150, reltol = 1e-6, tol = 1e-6)
#	
#EWMA.get.cc.Conditional.MC(p0 = 0.1, interval = c(2, 4), xmin = 0, xmax = 1, 
#	ymin = 0, ymax = 1, lambda = 0.2, eplison = 0, ARL0 = 370, mm = 50, ss = Inf, tt = 200, reltol = 1e-6, tol = 1e-6)

####################################################################################################################################################
#U <- 0.5
#V <- 0.5

#
#EWMA.CARL.MC.Q1 <- function(U, V, L, lambda, mm, ss, tt){
#	n <- 5
#	mu <- 0
#	sigma <- sqrt(n)
#	delta <- 0
#	#tt <- 100
#	z <- seq(-tt, tt, 1)
#	UCL <- L * sqrt(lambda / (2 - lambda))
#	LCL <- -UCL
#	v <- 2 * tt + 1
#	tau <- UCL / v
#	
#	S <- LCL + (2 * (tt + z) + 1) * tau
#	
#	Y <- matrix(0, 1, v)
#	Y[tt/2 + 1] <- 1
#	P <- diag(v)
#	I <- diag(v)
#	Q <- sqrt(qchisq(V, mm * (n - 1)) / mm / (n - 1))
#	Z <- qnorm(U)
#	#Q <- sqrt(rchisq(br, m * (n - 1)) / (m * (n - 1)))
#	#Z <- rnorm(br)
#	
#	for (l in 1:v){
#		#for (k in 1:v){
#			A <- Q * ((S + tau - (1 - lambda) * S[l]) / lambda) - ((sqrt(n) * delta / sigma) + (Z / sqrt(mm)))
#			B <- Q * ((S - tau - (1 - lambda) * S[l]) / lambda) - ((sqrt(n) * delta / sigma) + (Z / sqrt(mm)))
#			P[l, ] <- pnorm(A, mu, 1) - pnorm(B, mu, 1)
#		#}
#	}
#		
#	Y %*% solve(I - P) %*% matrix(1, v, 1)
#
#}
#

EWMA.CARL.MC.Q1 <- function(U, V, L, lambda, mm, ss, tt) {

	nn <- 2 * tt + 1

	Q <- matrix(NA, ncol = nn, nrow = nn)
	
	zz <- seq(-tt, tt, 1)
	
	for (l in 1:nn){
	
		ll <- l - (tt + 1)
	
		#for (k in 1:nn){
		
			#kk <- k - (tt + 1)
		
			#Q[l, k] <- pnorm(qnorm(U) / sqrt(mm) + (2 * kk - (1 - lambda) * 2 * ll + 1) / nn * 
			#	(L * sqrt(1 / lambda / (2 - lambda) * (1 - (1 - lambda)^(2 * ss))) * sqrt(qchisq(V, mm * (5 - 1))) / c4.f(mm * (5 - 1)) / sqrt(mm * (5 - 1)) )) - 
			#	pnorm(qnorm(U) / sqrt(mm) + (2 * kk - (1 - lambda) * 2 * ll - 1) / nn * 
			#	(L * sqrt(1 / lambda / (2 - lambda) * (1 - (1 - lambda)^(2 * ss))) * sqrt(qchisq(V, mm * (5 - 1))) / c4.f(mm * (5 - 1)) / sqrt(mm * (5 - 1)) )) 
				
				
			Q[l, ] <- pnorm(qnorm(U) / sqrt(mm) + (2 * zz - (1 - lambda) * 2 * ll + 1) / nn * 
				(L * sqrt(1 / lambda / (2 - lambda) * (1 - (1 - lambda)^(2 * ss))) * sqrt(qchisq(V, mm * (5 - 1))) / sqrt(mm * (5 - 1)) )) - 
				pnorm(qnorm(U) / sqrt(mm) + (2 * zz - (1 - lambda) * 2 * ll - 1) / nn * 
				(L * sqrt(1 / lambda / (2 - lambda) * (1 - (1 - lambda)^(2 * ss))) * sqrt(qchisq(V, mm * (5 - 1))) / sqrt(mm * (5 - 1)) )) 
		
		#}
	
	}
	
	I.matrix <- diag(nn)
	one.vec <- matrix(1, ncol = 1, nrow = nn)
	xi.vec <- matrix(0, ncol = nn, nrow = 1)
	xi.vec[1, tt / 2 + 1] <- 1
	
	xi.vec %*% solve(I.matrix - Q) %*% one.vec
	
}


EWMA.CARL.MC.Q2 <- function(U, V, L, lambda, mm, ss, tt){
	n <- 5
	mu <- 0
	sigma <- sqrt(n)
	delta <- 0
	#tt <- 100
	z <- seq(-tt, tt, 1)
	UCL <- L * sqrt(lambda / (2 - lambda))
	LCL <- -UCL
	v <- 2 * tt + 1
	tau <- UCL / v
	
	S <- LCL + (2 * (tt + z ) + 1) * tau
	#S <- numeric(length(z))
	#for (j in 1:length(z)){
	#	S[j] <- LCL + (2 * (v + z[j] + 1) * tau)
	#}
	
	Y <- matrix(0, 1, v)
	Y[tt/2 + 1] <- 1
	P <- diag(v)
	I <- diag(v)
	Q <- sqrt(qchisq(V, mm * (n - 1)) / mm / (n - 1))
	Z <- qnorm(U)
	
	for (l in 1:v){
		#for (k in 1:v){
			A <- Q * ((S + tau - (1 - lambda) * S[l]) / lambda) - ((sqrt(n) * delta / sigma) + (Z / sqrt(mm)))
			B <- Q * ((S - tau - (1 - lambda) * S[l]) / lambda) - ((sqrt(n) * delta / sigma) + (Z / sqrt(mm)))
			P[l, ] <- pnorm(A, mu, 1) - pnorm(B, mu, 1)
		#}
	}
		
	Y %*% solve(I - P) %*% matrix(1, v, 1)

	
}


###################################################################################################

#debug(EWMA.CARL.MC.Q1)
#EWMA.CARL.MC.Q1(U, V, 3.065, lambda = 0.2, mm = 300, ss = Inf, tt = 3)

#EWMA.CARL.MC.Q2(U, V, 3.065, lambda = 0.2, mm = 300, ss = Inf, tt = 3)


EWMA.CARL.Conditional.MC.integrand1 <- function(U, V, L, lambda, eplison, ARL0, mm, ss, tt){

	QQ <- EWMA.CARL.MC.Q1(U, V, L, lambda, mm, ss, tt)
	#xi.vec <- EWMA.CARL.MC.xi(tt)
	#xi.vec <- matrix(xi.vec, nrow = 1, ncol = 2 * tt + 1)
	#I.matrix <- diag(2 * tt + 1)
	#
	#one.vec <- matrix(1, ncol = 1, nrow = 2 * tt + 1)
	#
	#inv.matrix <- try(solve(I.matrix - QQ), silent = TRUE)
	#
	#if (class(inv.matrix) == 'try-error') {
	#
	#	cat('try-error', '\n')
	#	inv.matrix <- MASS::ginv(I.matrix - QQ)
	#
	#}
	#
	#ifelse(xi.vec %*% inv.matrix %*% one.vec >= (1 - eplison) * ARL0, 1, 0)
	
	ifelse(QQ >= (1 - eplison) * ARL0, 1, 0)

}

EWMA.CARL.Conditional.MC.integrand2 <- function(U, V, L, lambda, eplison, ARL0, mm, ss, tt){

	QQ <- EWMA.CARL.MC.Q2(U, V, L, lambda, mm, ss, tt)
	#xi.vec <- EWMA.CARL.MC.xi(tt)
	#xi.vec <- matrix(xi.vec, nrow = 1, ncol = 2 * tt + 1)
	#I.matrix <- diag(2 * tt + 1)
	#
	#one.vec <- matrix(1, ncol = 1, nrow = 2 * tt + 1)
	#
	#inv.matrix <- try(solve(I.matrix - QQ), silent = TRUE)
	#
	#if (class(inv.matrix) == 'try-error') {
	#
	#	cat('try-error', '\n')
	#	inv.matrix <- MASS::ginv(I.matrix - QQ)
	#
	#}
	#
	#ifelse(xi.vec %*% inv.matrix %*% one.vec >= (1 - eplison) * ARL0, 1, 0)
	
	ifelse(QQ >= (1 - eplison) * ARL0, 1, 0)

}

#debug(EWMA.CARL.MC)

#EWMA.CARL.MC(0.5, 0.5, 3, 0.2, 100, Inf, 5)

EWMA.CARL.Conditional.MC.integral1 <- function(xmin, xmax, ymin, ymax, L, lambda, eplison, ARL0, mm, ss, tt = 3, reltol = 1e-6){

	#cat(m, '\n')

	integral2(EWMA.CARL.Conditional.MC.integrand1, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, 
		singular = FALSE, vectorized = FALSE, L = L, lambda = lambda, eplison = eplison, ARL0 = ARL0, 
		mm = mm, ss = ss, tt = tt, reltol = reltol)$Q

}

EWMA.CARL.Conditional.MC.integral2 <- function(xmin, xmax, ymin, ymax, L, lambda, eplison, ARL0, mm, ss, tt = 3, reltol = 1e-6){

	#cat(m, '\n')

	integral2(EWMA.CARL.Conditional.MC.integrand2, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, 
		singular = FALSE, vectorized = FALSE, L = L, lambda = lambda, eplison = eplison, ARL0 = ARL0, 
		mm = mm, ss = ss, tt = tt, reltol = reltol)$Q

}



#debug(EWMA.CARL.MC.Q)

EWMA.get.cc.Conditional.MC1 <- function(p0 = 0.05, interval = c(1, 5), xmin = 0, xmax = 1, 
	ymin = 0, ymax = 1, lambda = 0.2, eplison = 0.1, ARL0 = 370, mm = 100, ss = Inf, tt = 3, reltol = 1e-6, tol = 1e-6){

		root.finding <- function(p0, xmin, xmax, ymin, ymax, L, lambda, eplison, ARL0, mm, ss, tt, reltol){
		
			p <- EWMA.CARL.Conditional.MC.integral1(xmin, xmax, ymin, ymax, 
				L, lambda, eplison, ARL0, mm, ss, tt, reltol)
			
			cat('L:', L, ', p:', 1 - p, '\n')
		
			1 - p0 - p
		
		}
	
		uniroot(root.finding, interval = interval, tol = tol, p0 = p0, xmin = xmin, xmax = xmax,
			ymin = ymin, ymax = ymax, lambda = lambda, eplison = eplison, ARL0 = ARL0, 
			mm = mm, ss = ss, tt = tt, reltol = reltol)$root
	
}

EWMA.get.cc.Conditional.MC2 <- function(p0 = 0.05, interval = c(1, 5), xmin = 0, xmax = 1, 
	ymin = 0, ymax = 1, lambda = 0.2, eplison = 0.1, ARL0 = 370, mm = 100, ss = Inf, tt = 3, reltol = 1e-6, tol = 1e-6){

		root.finding <- function(p0, xmin, xmax, ymin, ymax, L, lambda, eplison, ARL0, mm, ss, tt, reltol){
		
			p <- EWMA.CARL.Conditional.MC.integral2(xmin, xmax, ymin, ymax, 
				L, lambda, eplison, ARL0, mm, ss, tt, reltol)
			
			cat('L:', L, ', p:', 1 - p, '\n')
		
			1 - p0 - p
		
		}
	
		uniroot(root.finding, interval = interval, tol = tol, p0 = p0, xmin = xmin, xmax = xmax,
			ymin = ymin, ymax = ymax, lambda = lambda, eplison = eplison, ARL0 = ARL0, 
			mm = mm, ss = ss, tt = tt, reltol = reltol)$root
	
}


#EWMA.get.cc.Conditional.MC1(p0 = 0.1, interval = c(2.3, 4), xmin = 0, xmax = 1, 
#	ymin = 0, ymax = 1, lambda = 0.2, eplison = 0.1, ARL0 = 370, mm = 300, ss = Inf, tt = 50, reltol = 1e-6, tol = 1e-6)
#	
#EWMA.get.cc.Conditional.MC1(p0 = 0.1, interval = c(2.3, 4), xmin = 0, xmax = 1, 
#	ymin = 0, ymax = 1, lambda = 0.2, eplison = 0.1, ARL0 = 370, mm = 300, ss = Inf, tt = 100, reltol = 1e-6, tol = 1e-6)
#	
#EWMA.get.cc.Conditional.MC1(p0 = 0.1, interval = c(2.3, 4), xmin = 0, xmax = 1, 
#	ymin = 0, ymax = 1, lambda = 0.2, eplison = 0.1, ARL0 = 370, mm = 300, ss = Inf, tt = 150, reltol = 1e-6, tol = 1e-6)
#	
#EWMA.get.cc.Conditional.MC1(p0 = 0.1, interval = c(2.3, 4), xmin = 0, xmax = 1, 
#	ymin = 0, ymax = 1, lambda = 0.2, eplison = 0.1, ARL0 = 370, mm = 300, ss = Inf, tt = 200, reltol = 1e-6, tol = 1e-6)

#EWMA.get.cc.Conditional.MC2(p0 = 0.1, interval = c(2.3, 4), xmin = 0, xmax = 1, 
#	ymin = 0, ymax = 1, lambda = 0.2, eplison = 0.1, ARL0 = 370, mm = 300, ss = Inf, tt = 20, reltol = 1e-6, tol = 1e-6)
#
