effectiveAlpha.mvNorm <- function(rho, k, n = 11) {
	if(length(rho) > 1) {
		out <- numeric(length(rho))
		for(i in 1:length(rho)) {
			out[i] <- effectiveAlpha.mvNorm(rho[i], k, n)
		}
	} else if(length(k) > 1) {
		out <- numeric(length(k))
		for(i in 1:length(k)) {
			out[i] <- effectiveAlpha.mvNorm(rho, k[i], n)
		}
	} else {
		sigma <- diag(rep(1, n))
		sigma <- rho ^ abs(row(sigma) - col(sigma))
		out <- 1 - pmvnorm(-k, k, mean = rep(0, n), sigma)[1]
	}
	
	out
}