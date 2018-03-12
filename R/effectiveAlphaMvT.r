effectiveAlpha.mvT <- function(rho, k, n = 11, df = 10) {
	if(length(rho) > 1) {
		out <- numeric(length(rho))
		for(i in 1:length(rho)) {
			out[i] <- effectiveAlpha.mvT(rho[i], k, n, df)
		}
	} else if(length(k) > 1) {
		out <- numeric(length(k))
		for(i in 1:length(k)) {
			out[i] <- effectiveAlpha.mvT(rho, k[i], n, df)
		}
	} else {
		sigma <- diag(rep(1, n))
		sigma <- rho ^ abs(row(sigma) - col(sigma))
		out <- 1 - pmvt(-k, k, delta = rep(0, n), df = df, sigma)[1]
	}
	
	out
}