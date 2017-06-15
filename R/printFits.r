printFits <- function(part1.list) {
	R2.g1.1 <- part1.list$R2.g1.1
	R2.g1.2 <- part1.list$R2.g1.2
	R2.g2.1 <- part1.list$R2.g2.1
	R2.g2.2 <- part1.list$R2.g2.2
	
	cor.1 <- part1.list$cor.1
	cor.2 <- part1.list$cor.2
	cor.3 <- part1.list$cor.3
	cor.4 <- part1.list$cor.4
	
	diffs <- part1.list$diffs

	ar1.good <- sum(c(R2.g1.1, R2.g1.2, R2.g2.1, R2.g2.2) >= .95 &
		c(cor.1, cor.3, cor.2, cor.4), na.rm = TRUE)
	ar1.ok <- sum(c(R2.g1.1, R2.g1.2, R2.g2.1, R2.g2.2) < .95 &
		c(R2.g1.1, R2.g1.2, R2.g2.1, R2.g2.2) >= .8 &
		c(cor.1, cor.3, cor.2, cor.4), na.rm = TRUE)
	ar1.bad <- sum(c(R2.g1.1, R2.g1.2, R2.g2.1, R2.g2.2) < .8 &
		c(R2.g1.1, R2.g1.2, R2.g2.1, R2.g2.2) > 0 &
		c(cor.1, cor.3, cor.2, cor.4), na.rm = TRUE)
		
	nonar1.good <- sum(c(R2.g1.1, R2.g1.2, R2.g2.1, R2.g2.2) >= .95 &
		!c(cor.1, cor.3, cor.2, cor.4), na.rm = TRUE)
	nonar1.ok <- sum(c(R2.g1.1, R2.g1.2, R2.g2.1, R2.g2.2) < .95 &
		c(R2.g1.1, R2.g1.2, R2.g2.1, R2.g2.2) >= .8 &
		!c(cor.1, cor.3, cor.2, cor.4), na.rm = TRUE)
	nonar1.bad <- sum(c(R2.g1.1, R2.g1.2, R2.g2.1, R2.g2.2) < .8 &
		c(R2.g1.1, R2.g1.2, R2.g2.1, R2.g2.2) > 0 &
		!c(cor.1, cor.3, cor.2, cor.4), na.rm = TRUE)
		
	if(diffs) {
		nofit <- sum(is.na(c(cor.1, cor.3, cor.2, cor.4)))
	} else {
		nofit <- sum(is.na(c(cor.1, cor.2)))
	}

	cat("########################################\n")
	cat("############### FITS ###################\n")
	cat("########################################\n")
	cat(paste("AR1,       R2>=0.95   --", ar1.good, "\n"))
	cat(paste("AR1,     0.95>R2>=0.8 --", ar1.ok, "\n"))
	cat(paste("AR1,       0.8>R2     --", ar1.bad, "\n"))
	cat(paste("Non-AR1,   R2>=0.95   --", nonar1.good, "\n"))
	cat(paste("Non-AR1, 0.95>R2>=0.8 --", nonar1.ok, "\n"))
	cat(paste("Non-AR1,   0.8>R2     --", nonar1.bad, "\n"))
	cat(paste("No Fit                --", nofit, "\n"))
	cat("########################################\n\n")
}
