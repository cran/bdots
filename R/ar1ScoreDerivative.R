ar1.score.derivative <- function(phi, y) {
	val.1 <- (1 + phi ^ 2) / (1 - phi ^ 2) ^ 2 * length(y)
	val.2 <- (1 + 3 * phi ^ 2) / (1 - phi ^ 2) ^ 3 * (2 * sum(y[1:(length(y) - 1)] ^ 2) + y[length(y)] ^ 2)
	val.3 <- 2 * phi * (3 + phi ^ 2) / (1 - phi ^ 2) ^ 3 * sum(y[2:length(y)] * y[1:(length(y) - 1)])
	
	val.1 - val.2 + val.3
}