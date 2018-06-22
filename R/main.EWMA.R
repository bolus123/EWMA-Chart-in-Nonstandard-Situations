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
	denom <- lambda / (2 - lambda) / (1 - (1 - lambda) ^ (2 * ss))
	
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

EWMA.get.cc <- function(ARL0 = 370, interval = c(1, 5), xmin = 0, xmax = 1, 
	ymin = 0, ymax = 1, lambda = 0.2, mm = 100, ss = Inf, reltol = 1e-6){

		root.finding <- function(ARL0, xmin, xmax, ymin, ymax, L, lambda, mm, ss, reltol){
		
			ARL0 - EWMA.CARL.DP.integral(xmin, xmax, ymin, ymax, L, lambda, mm, ss, reltol)
		
		}
	
		uniroot(root.finding, interval = interval, ARL0 = ARL0, xmin = xmin, xmax = xmax,
			ymin = ymin, ymax = ymax, lambda = lambda, mm = mm, ss = ss, reltol = reltol)$root
	
}


aa <- EWMA.get.cc(ARL0 = 370, interval = c(1, 2.57), xmin = 0, xmax = 1, 
	ymin = 0, ymax = 1, lambda = 0.2, mm = 10, ss = Inf)
	
bb <- EWMA.get.cc(ARL0 = 370, interval = c(1, 3), xmin = 0, xmax = 1, 
	ymin = 0, ymax = 1, lambda = 0.2, mm = 20, ss = Inf)
	
cc <- EWMA.get.cc(ARL0 = 370, interval = c(1, 4), xmin = 0, xmax = 1, 
	ymin = 0, ymax = 1, lambda = 0.2, mm = 50, ss = Inf)
	
dd <- EWMA.get.cc(ARL0 = 370, interval = c(1, 4), xmin = 0, xmax = 1, 
	ymin = 0, ymax = 1, lambda = 0.2, mm = 100, ss = Inf)	
	
ee <- EWMA.get.cc(ARL0 = 370, interval = c(1, 4), xmin = 0, xmax = 1, 
	ymin = 0, ymax = 1, lambda = 0.2, mm = 200, ss = Inf)	
	
ff <- EWMA.get.cc(ARL0 = 370, interval = c(1, 4), xmin = 0, xmax = 1, 
	ymin = 0, ymax = 1, lambda = 0.2, mm = 300, ss = Inf)	

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

EWMA.ARLin.performance(L = 2.963359, lambda = 0.2, mm = 20, ss = Inf, sim = 1000000)

####################################################################################################################################################
