ar1.solve  <- function(y) {
	val.1 <- - length(y)
	val.2 <- sum(y[2:length(y)] * y[1:(length(y) - 1)])
	val.3 <- length(y) - 2 * sum(y[1:(length(y) - 1)] ^ 2) - y[length(y)] ^ 2
	val.4 <- val.2
	
	phi.est <- polyroot(c(val.4, val.3, val.2, val.1))
	phi.est <- Re(subset(phi.est, abs(Im(phi.est)) < 0.0001)) # Discard complex roots
	phi.est <- subset(phi.est, phi.est >= 0)
	phi.est <- subset(phi.est, ar1.score.derivative(phi.est, y) < 0) # Discard roots where our function is not concave
	if(length(phi.est) > 1) phi.est <- phi.est[which.max(ar1.loglikelihood(phi.est, y))]
	if(length(phi.est) == 0) {
		warning("AR1 solver could not find a proper root. Using Yule-Walker solver.")
		phi.est <- ar(y, order.max = 1, aic = FALSE)$ar
	}

	phi.est
}