ar1.loglikelihood <- function(phi, y) {
	n <- length(y)
	val.1 <- - n / 2 * log(2 * pi * (1 - phi ^ 2))
	val.2 <- y[1] ^ 2 + sum((y[2:n] - phi * y[1:(n - 1)]) ^ 2)
	val.3 <- 2 * (1 - phi ^ 2)
	
	val.1 - val.2 / val.3
}