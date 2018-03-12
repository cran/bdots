effectiveAlpha.normApprox <- function(rho, k, n = 11) {
	if(length(rho) > 1) {
		out <- numeric(length(rho))
		for(i in 1:length(rho)) {
			out[i] <- effectiveAlpha.normApprox(rho[i], k, n)
		}
	} else if(length(k) > 1) {
		out <- numeric(length(k))
		for(i in 1:length(k)) {
			out[i] <- effectiveAlpha.normApprox(rho, k[i], n)
		}
	} else {
		p1 <- 2 * pnorm(k) - 1
		mean <- c(0, 0)
		sigma <- matrix(c(1, rho, rho, 1), ncol = 2)
		p2 <- pmvnorm(c(-k, -k), c(k, k), mean = mean, sigma = sigma, algorithm=GenzBretz(abseps=.Machine$double.eps))[1]
		p3 <- p2 / p1
		out <- 1 - p3 ^ (n - 1) * p1
	}
	
	out
}