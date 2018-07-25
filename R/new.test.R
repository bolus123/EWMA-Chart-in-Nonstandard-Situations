library(mvtnorm)

variance.matrix.f <- function(lambda, sigma2, n){

	first.row <- unlist(lapply(
		1:n,
		function(x) lambda / (2 - lambda) * sigma2 * (1 - lambda) ^ (x - 1)
	))

	toeplitz(first.row)

}

run.length.prob.f <- function(n, cc, mu, sigma2, lambda) {

	mu.vec <- rep(mu, n)
	var.mat <- variance.matrix.f(lambda, sigma2, n)
	
	lower1 <- c(rep(-cc, n - 1), -Inf)
	upper1 <- c(rep(cc, n - 1), -cc)
	
	lower2 <- c(rep(-cc, n - 1), cc)
	upper2 <- c(rep(cc, n - 1), Inf)
	
	pmvnorm(lower = lower1, upper = upper1, mean=mu.vec, sigma=var.mat) + 
	pmvnorm(lower = lower2, upper = upper2, mean=mu.vec, sigma=var.mat)

}

a <- rep(NA, 50)

for (i in 1:50) {

	a[i] <- run.length.prob.f(i, 0.5, 0, 1, 0.2)

}