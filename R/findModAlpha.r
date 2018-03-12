find.mod.alpha <- function(rho, n, df, alpha = .05, error.acc = .001, grad.diff = ifelse(cores > 3, .5, .1), verbose = FALSE, cores = 2, method = "t") {
	if(n > 1000) stop("Too many testing locations to accurately calculate modified alpha")
	if(method != "NormApprox") f <- function(k) effectiveAlpha.mvT(rho, k, n, df) else f <- function(k) effectiveAlpha.mvNorm(rho, k, n)
	if(cores > 1 && method != "old") {
		cl <- makePSOCKcluster(rep("localhost", cores))
		clusterExport(cl, c("f", "effectiveAlpha.mvT", "effectiveAlpha.mvNorm", "rho", "n", "df", "pmvt", "pmvnorm"), envir = environment())
	}
	
	min.val <- qt(1 - alpha / 2, df)
	max.val <- qt(1 - alpha / 2, df) * 2
	while(effectiveAlpha.normApprox(rho, max.val, n) >= alpha) max.val <- max.val * 2
	k.f1 <- uniroot(function(k) effectiveAlpha.normApprox(rho, k, n) - alpha, interval = c(min.val, max.val))$root

	k.f2 <- k.f1
	if (method != "old") { # Second - use Normal Approximation estimate to find the root using Multivariate Normal, can use slope approximation of log response
		
		if(cores > 3) {
			if(cores %% 2 == 1) {
				spots <- seq(k.f2 - grad.diff, k.f2 + grad.diff, length.out = cores)
				spots <- c(k.f2, spots[-(cores + 1) / 2])
			} else {
				spots <- seq(k.f2 - grad.diff, k.f2 + grad.diff, length.out = cores)
				spots <- c(k.f2, spots[-sample(cores, 1)])
			}
			alpha.star.full <- parLapply(cl, spots, function(k) f(k))
			alpha.star <- alpha.star.full[[1]]
		} else if(cores == 3) {
			alpha.star.full <- parLapply(cl, c(k.f2, k.f2 + grad.diff, k.f2 - grad.diff), function(k) f(k))
			alpha.star <- alpha.star.full[[1]]
			alpha.star.2 <- alpha.star.full[[2]]
			alpha.star.3 <- alpha.star.full[[3]]
		}else if(cores == 2) {
			alpha.star.full <- parLapply(cl, c(k.f2, k.f2 + grad.diff), function(k) f(k))
			alpha.star <- alpha.star.full[[1]]
			alpha.star.2 <- alpha.star.full[[2]]
		} else {
			alpha.star <- f(k.f2)
			alpha.star.2 <- f(k.f2 + grad.diff)
			alpha.star.full <- list(alpha.star, alpha.star.2)
		}
		if(verbose) print(paste("Bivariate approximation for alpha:", alpha.star, ", for k:", k.f2))

		error.min <- min(abs(unlist(alpha.star.full) - alpha))
		if(error.min <= error.acc) alpha.star <- alpha.star.full[[which(abs(unlist(alpha.star.full) - alpha) == error.min)]]
		
		while(abs(alpha.star - alpha) > error.acc) {
			if(cores > 3){ # Fit quadratic linear model and use quadratic formula to find lower root
				spots2 <- spots ^ 2
				fit <- lm(unlist(alpha.star.full) ~ spots + spots2)
				coef <- coef(fit)
				print(coef)
				k.f2 <- as.vector((-coef[2] - sqrt(coef[2] ^ 2 - 4 * coef[3] * (coef[1] - alpha))) / (2 * coef[3]))
				
				if(cores %% 2 == 1) {
					spots <- seq(k.f2 - grad.diff, k.f2 + grad.diff, length.out = cores)
					spots <- c(k.f2, spots[-(cores + 1) / 2])
				} else {
					spots <- seq(k.f2 - grad.diff, k.f2 + grad.diff, length.out = cores)
					spots <- c(k.f2, spots[-sample(cores, 1)])
				}
				alpha.star.full <- parLapply(cl, spots, function(k) f(k))
				alpha.star <- alpha.star.full[[1]]
			} else if(cores == 3) { # Use Halley's Method (i.e. estimate first two derivatives) -- CHECK THIS FORMULA, NOT WORKING CORRECTLY
				grad.est <- (log(alpha.star.2) - log(alpha.star.3)) / (2 * grad.diff)
				hess.est <- (log(alpha.star.2) - 2 * log(alpha.star) + log(alpha.star.3)) / grad.diff ^ 2
				k.f2 <- k.f2 - (log(alpha.star) - log(alpha)) / grad.est * (1 - (log(alpha.star) - log(alpha)) * hess.est / (2 * grad.est ^ 2)) ^ (-1)
				
				alpha.star.full <- parLapply(cl, c(k.f2, k.f2 + grad.diff, k.f2 - grad.diff), function(k) f(k))
				alpha.star <- alpha.star.full[[1]]
				alpha.star.2 <- alpha.star.full[[2]]
				alpha.star.3 <- alpha.star.full[[3]]
			} else { # Use Newton's Method (i.e. estimate gradient and find root from that)
				grad.est <- (log(alpha.star.2) - log(alpha.star)) / grad.diff
				k.f2 <- k.f2 - (log(alpha.star) - log(alpha)) / grad.est
				
				if(cores == 2) { 
					alpha.star.full <- parLapply(cl, c(k.f2, k.f2 + grad.diff), function(k) f(k))
					alpha.star <- alpha.star.full[[1]]
					alpha.star.2 <- alpha.star.full[[2]]
				} else {
					alpha.star <- f(k.f2)
					alpha.star.2 <- f(k.f2 + grad.diff)
					alpha.star.full <- list(alpha.star, alpha.star.2)
				}
			}
			
			if(verbose) print(paste("Full approximation for alpha:", alpha.star, ", for k:", k.f2))
			
			error.min <- min(abs(unlist(alpha.star.full) - alpha))
			if(error.min <= error.acc) alpha.star <- alpha.star.full[[which(abs(unlist(alpha.star.full) - alpha) == error.min)]]
		}
		
		if(cores > 1) stopCluster(cl)
	}

	out <- pt(k.f2, df, lower.tail = FALSE) * 2
  
  out
}