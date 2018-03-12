doubleGauss.boot <- function(part1.list, seed = new.seed(), alpha = 0.05, paired = FALSE, N.iter = 1000, cores = 1, p.adj = "oleson", test.spots = NULL, time.test = NULL, test.params = FALSE) {
  data       <- part1.list$data
  col        <- part1.list$col
  rho.0      <- part1.list$rho.0
  N.time     <- part1.list$N.time
  coef.id1   <- part1.list$coef.id1
  coef.id2   <- part1.list$coef.id2
	coef.id3   <- part1.list$coef.id3
  coef.id4   <- part1.list$coef.id4
	sdev.id1.cov   <- part1.list$sdev.id1.cov
	sdev.id2.cov   <- part1.list$sdev.id2.cov
	sdev.id3.cov   <- part1.list$sdev.id3.cov
	sdev.id4.cov   <- part1.list$sdev.id4.cov
  sigma.id1  <- part1.list$sigma.id1
  sigma.id2  <- part1.list$sigma.id2
	sigma.id3  <- part1.list$sigma.id3
  sigma.id4  <- part1.list$sigma.id4
  id.nums.g1 <- part1.list$id.nums.g1
  id.nums.g2 <- part1.list$id.nums.g2
  groups     <- part1.list$groups
  time.all   <- part1.list$time.all
  N.g1       <- part1.list$N.g1
  N.g2       <- part1.list$N.g2
	diffs      <- part1.list$diffs
	
	if(!is.null(test.spots)) time.all <- test.spots
	N.tests <- length(time.all)
	N.time <- length(time.all)
	
  if(diffs) {
    group1.bad <- is.na(coef.id1[,1]) | is.na(coef.id3[,1])
    group2.bad <- is.na(coef.id2[,1]) | is.na(coef.id4[,1])
  } else {
    group1.bad <- is.na(coef.id1[,1])
    group2.bad <- is.na(coef.id2[,1])
  }
	
	coef.id1 <- subset(coef.id1, !group1.bad)
	coef.id3 <- subset(coef.id3, !group1.bad)
	coef.id2 <- subset(coef.id2, !group2.bad)
	coef.id4 <- subset(coef.id4, !group2.bad)
	
	sdev.id1.cov <- subset(sdev.id1.cov, !group1.bad)
	if(diffs) sdev.id3.cov <- subset(sdev.id3.cov, !group1.bad)
	sdev.id2.cov <- subset(sdev.id2.cov, !group2.bad)
	if(diffs) sdev.id4.cov <- subset(sdev.id4.cov, !group2.bad)
	
	sigma.id1 <- subset(sigma.id1, !group1.bad)
	sigma.id3 <- subset(sigma.id3, !group1.bad)
	sigma.id2 <- subset(sigma.id2, !group2.bad)
	sigma.id4 <- subset(sigma.id4, !group2.bad)
	
	id.nums.g1 <- id.nums.g1[!group1.bad]
	id.nums.g2 <- id.nums.g2[!group2.bad]
	N.g1 <- N.g1 - sum(group1.bad)
	N.g2 <- N.g2 - sum(group2.bad)

	curve1.0 <-  matrix(NA, ncol = N.time, nrow = N.g1)
	mu1.ran <-   rep(NA, N.g1)
	ht1.ran <-   rep(NA, N.g1)
	s11.ran <-   rep(NA, N.g1)
	s21.ran <-   rep(NA, N.g1)
	b11.ran <-   rep(NA, N.g1)
	b21.ran <-   rep(NA, N.g1)

	curve2.0 <-  matrix(NA, ncol = N.time, nrow = N.g2)
	mu2.ran <-   rep(NA, N.g2)
	ht2.ran <-   rep(NA, N.g2)
	s12.ran <-   rep(NA, N.g2)
	s22.ran <-   rep(NA, N.g2)
	b12.ran <-   rep(NA, N.g2)
	b22.ran <-   rep(NA, N.g2)

	curve1.mat <- matrix(NA, ncol = N.time, nrow = N.iter)
	curve2.mat <- matrix(NA, ncol = N.time, nrow = N.iter)
	curve3.mat <- matrix(NA, ncol = N.time, nrow = N.iter)
	
	if(paired && N.g1 != N.g2) stop("Can't do paired test for different group sizes")
  
  #Target Curve
	curve.f <- function(mu, ht, sig1, sig2, base1, base2, x){
	  whichgauss <- x < mu
	  y1 <- exp(-1 * (x - mu) ^ 2 / (2 * sig1 ^ 2)) * (ht - base1) + base1
	  y2 <- exp(-1 * (x - mu) ^ 2 / (2 * sig2 ^ 2)) * (ht - base2) + base2
	
	  y <- whichgauss * y1 + (1 - whichgauss) * y2
		y
	}
	
	mu.1 <- ht.1 <- s1.1 <- s2.1 <- b1.1 <- b2.1 <-
		mu.2 <- ht.2 <- s1.2 <- s2.2 <- b1.2 <- b2.2 <- numeric(N.iter)
		
	# Generate parameters
	params.1 <- vector("list", N.g1)
	params.2 <- vector("list", N.g2)
	if(diffs) {
		params.3 <- params.1
		params.4 <- params.2
	}
	
	if(paired) {
		param.cor <- cor(coef.id1, coef.id2)
		
		for(i in 1:N.g1) {
			sigma.11 <- sdev.id1.cov[[i]]
			sigma.22 <- sdev.id2.cov[[i]]
			sigma.12 <- param.cor * sqrt(t(t(diag(sigma.11))) %*% t(diag(sigma.22)))
			sigma <- cbind(rbind(sigma.11, t(sigma.12)), rbind(sigma.12, sigma.22))
			while(any(eigen(sigma)$values <= 0)) sigma <- as.matrix(nearPD(sigma, keepDiag = TRUE, maxit = 1e7)$mat)
			
			param.i <- rmvnorm(N.iter, c(coef.id1[i,], coef.id2[i,]), sigma)
			params.1[[i]] <- param.i[,1:6]
			params.2[[i]] <- param.i[,7:12]
		}
		
		if(diffs) {
			param.cor <- cor(coef.id3, coef.id4)
			
			for(i in 1:N.g1) {
				sigma.11 <- sdev.id3.cov[[i]]
				sigma.22 <- sdev.id4.cov[[i]]
				sigma.12 <- param.cor * sqrt(t(t(diag(sigma.11))) %*% t(diag(sigma.22)))
				sigma <- cbind(rbind(sigma.11, t(sigma.12)), rbind(sigma.12, sigma.22))
				while(any(eigen(sigma)$values <= 0)) sigma <- as.matrix(nearPD(sigma, keepDiag = TRUE, maxit = 1e7)$mat)
				
				param.i <- rmvnorm(N.iter, c(coef.id1[i,], coef.id2[i,]), sigma)
				params.3[[i]] <- param.i[,1:6]
				params.4[[i]] <- param.i[,7:12]
			}
		}
	} else {
		for(i in 1:N.g1) {
			params.1[[i]] <- rmvnorm(N.iter, coef.id1[i,], sdev.id1.cov[[i]])
			if(diffs) params.3[[i]] <- rmvnorm(N.iter, coef.id3[i,], sdev.id3.cov[[i]])
		}
		for(i in 1:N.g2) {
			params.2[[i]] <- rmvnorm(N.iter, coef.id2[i,], sdev.id2.cov[[i]])
			if(diffs) params.4[[i]] <- rmvnorm(N.iter, coef.id4[i,], sdev.id4.cov[[i]])
		}
	}
	
	##################
	##### 1 Core #####
	##################
	if(cores == 1) {
		set.seed(seed)
		for(iter in 1:N.iter){
				mu.1[iter] <- mean(vapply(params.1, function(x) x[iter, 1], numeric(1)))
				ht.1[iter] <- mean(vapply(params.1, function(x) x[iter, 2], numeric(1)))
				s1.1[iter] <- mean(vapply(params.1, function(x) x[iter, 3], numeric(1)))
				s2.1[iter] <- mean(vapply(params.1, function(x) x[iter, 4], numeric(1)))
				b1.1[iter] <- mean(vapply(params.1, function(x) x[iter, 5], numeric(1)))
				b2.1[iter] <- mean(vapply(params.1, function(x) x[iter, 6], numeric(1)))

				mu.2[iter] <- mean(vapply(params.1, function(x) x[iter, 1], numeric(1)))
				ht.2[iter] <- mean(vapply(params.1, function(x) x[iter, 2], numeric(1)))
				s1.2[iter] <- mean(vapply(params.1, function(x) x[iter, 3], numeric(1)))
				s2.2[iter] <- mean(vapply(params.1, function(x) x[iter, 4], numeric(1)))
				b1.2[iter] <- mean(vapply(params.1, function(x) x[iter, 5], numeric(1)))
				b2.2[iter] <- mean(vapply(params.1, function(x) x[iter, 6], numeric(1)))
				
			if(diffs) {
				for(id in 1:N.g1){ #Get fixation level for each g1 subject
					curve1.0[id,] <- curve.f(params.1[[id]][iter, 1],
																		params.1[[id]][iter, 2],
																		params.1[[id]][iter, 3],
																		params.1[[id]][iter, 4],
																		params.1[[id]][iter, 5],
																		params.1[[id]][iter, 6],
																		time.all) -
													 curve.f(params.3[[id]][iter, 1],
																		params.3[[id]][iter, 2],
																		params.3[[id]][iter, 3],
																		params.3[[id]][iter, 4],
																		params.3[[id]][iter, 5],
																		params.3[[id]][iter, 6],
																		time.all)
				}
				for(id in 1:N.g2){ #Get fixation level for each g2 subject
					curve2.0[id,] <- curve.f(params.2[[id]][iter, 1],
																		params.2[[id]][iter, 2],
																		params.2[[id]][iter, 3],
																		params.2[[id]][iter, 4],
																		params.2[[id]][iter, 5],
																		params.2[[id]][iter, 6],
																		time.all) -
													 curve.f(params.4[[id]][iter, 1],
																		params.4[[id]][iter, 2],
																		params.4[[id]][iter, 3],
																		params.4[[id]][iter, 4],
																		params.4[[id]][iter, 5],
																		params.4[[id]][iter, 6],
																		time.all)
				}
			} else {
				for(id in 1:N.g1){ #Get fixation level for each g1 subject
					curve1.0[id,] <- curve.f(params.1[[id]][iter, 1],
																		params.1[[id]][iter, 2],
																		params.1[[id]][iter, 3],
																		params.1[[id]][iter, 4],
																		params.1[[id]][iter, 5],
																		params.1[[id]][iter, 6],
																		time.all)
				}
				for(id in 1:N.g2){ #Get fixation level for each g2 subject
					curve2.0[id,] <- curve.f(params.2[[id]][iter, 1],
																		params.2[[id]][iter, 2],
																		params.2[[id]][iter, 3],
																		params.2[[id]][iter, 4],
																		params.2[[id]][iter, 5],
																		params.2[[id]][iter, 6],
																		time.all)
				}
			}

			curve1.mat[iter,] <- colMeans(curve1.0)
			curve2.mat[iter,] <- colMeans(curve2.0)
			
			if(paired) curve3.mat[iter,] <- curve2.mat[iter,] - curve1.mat[iter,]
		}
	} else {
	####################
	##### 2+ Cores #####
	####################
	
		cl <- makePSOCKcluster(cores)
		registerDoParallel(cl)
		
		for.out <- foreach(iter = 1:N.iter, .combine = rbind, .options.RNG = seed) %dorng% {
			mu.temp.1 <- mean(vapply(params.1, function(x) x[iter, 1], numeric(1)))
			ht.temp.1 <- mean(vapply(params.1, function(x) x[iter, 2], numeric(1)))
			s1.temp.1 <- mean(vapply(params.1, function(x) x[iter, 3], numeric(1)))
			s2.temp.1 <- mean(vapply(params.1, function(x) x[iter, 4], numeric(1)))
			b1.temp.1 <- mean(vapply(params.1, function(x) x[iter, 5], numeric(1)))
			b2.temp.1 <- mean(vapply(params.1, function(x) x[iter, 6], numeric(1)))

			mu.temp.2 <- mean(vapply(params.1, function(x) x[iter, 1], numeric(1)))
			ht.temp.2 <- mean(vapply(params.1, function(x) x[iter, 2], numeric(1)))
			s1.temp.2 <- mean(vapply(params.1, function(x) x[iter, 3], numeric(1)))
			s2.temp.2 <- mean(vapply(params.1, function(x) x[iter, 4], numeric(1)))
			b1.temp.2 <- mean(vapply(params.1, function(x) x[iter, 5], numeric(1)))
			b2.temp.2 <- mean(vapply(params.1, function(x) x[iter, 6], numeric(1)))
			
			if(diffs) {
				for(id in 1:N.g1){ #Get fixation level for each g1 subject
					curve1.0[id,] <- curve.f(params.1[[id]][iter, 1],
																		params.1[[id]][iter, 2],
																		params.1[[id]][iter, 3],
																		params.1[[id]][iter, 4],
																		params.1[[id]][iter, 5],
																		params.1[[id]][iter, 6],
																		time.all) -
													 curve.f(params.3[[id]][iter, 1],
																		params.3[[id]][iter, 2],
																		params.3[[id]][iter, 3],
																		params.3[[id]][iter, 4],
																		params.3[[id]][iter, 5],
																		params.3[[id]][iter, 6],
																		time.all)
				}
				for(id in 1:N.g2){ #Get fixation level for each g2 subject
					curve2.0[id,] <- curve.f(params.2[[id]][iter, 1],
																		params.2[[id]][iter, 2],
																		params.2[[id]][iter, 3],
																		params.2[[id]][iter, 4],
																		params.2[[id]][iter, 5],
																		params.2[[id]][iter, 6],
																		time.all) -
													 curve.f(params.4[[id]][iter, 1],
																		params.4[[id]][iter, 2],
																		params.4[[id]][iter, 3],
																		params.4[[id]][iter, 4],
																		params.4[[id]][iter, 5],
																		params.4[[id]][iter, 6],
																		time.all)
				}
			} else {
				for(id in 1:N.g1){ #Get fixation level for each g1 subject
					curve1.0[id,] <- curve.f(params.1[[id]][iter, 1],
																		params.1[[id]][iter, 2],
																		params.1[[id]][iter, 3],
																		params.1[[id]][iter, 4],
																		params.1[[id]][iter, 5],
																		params.1[[id]][iter, 6],
																		time.all)
				}
				for(id in 1:N.g2){ #Get fixation level for each g2 subject
					curve2.0[id,] <- curve.f(params.2[[id]][iter, 1],
																		params.2[[id]][iter, 2],
																		params.2[[id]][iter, 3],
																		params.2[[id]][iter, 4],
																		params.2[[id]][iter, 5],
																		params.2[[id]][iter, 6],
																		time.all)
				}
			}

			curve1 <- colMeans(curve1.0)
			curve2 <- colMeans(curve2.0)
			curve3 <- curve2 - curve1
			
			c(curve1, curve2, curve3, mu.temp.1, ht.temp.1, s1.temp.1, s2.temp.1, b1.temp.1, b2.temp.1,
				mu.temp.2, ht.temp.2, s1.temp.2, s2.temp.2, b1.temp.2, b2.temp.2)
		}
		
		curve1.mat <- for.out[,1:N.time]
		curve2.mat <- for.out[,(N.time + 1):(2 * N.time)]
		curve3.mat <- for.out[,(2 * N.time + 1):(3 * N.time)]
		mu.1 <- for.out[, 3 * N.time + 1]
		ht.1 <- for.out[, 3 * N.time + 2]
		s1.1 <- for.out[, 3 * N.time + 3]
		s2.1 <- for.out[, 3 * N.time + 4]
		b1.1 <- for.out[, 3 * N.time + 5]
		b2.1 <- for.out[, 3 * N.time + 6]
		mu.2 <- for.out[, 3 * N.time + 7]
		ht.2 <- for.out[, 3 * N.time + 8]
		s1.2 <- for.out[, 3 * N.time + 9]
		s2.2 <- for.out[, 3 * N.time + 10]
		b1.2 <- for.out[, 3 * N.time + 11]
		b2.2 <- for.out[, 3 * N.time + 12]
		
		stopCluster(cl)
	}
	
	curve.mean1 <- apply(curve1.mat, 2, mean)
	curve.mean2 <- apply(curve2.mat, 2, mean)
	
	curve.g1    <- curve.mean1
	curve.g2    <- curve.mean2
	
	curve.sd1   <- apply(curve1.mat, 2, sd)
	curve.sd2   <- apply(curve2.mat, 2, sd)

	if(paired) {
		diff.mean <- apply(curve3.mat, 2, mean)
		curve.sd <- apply(curve3.mat, 2, sd)
		
		t.val <- diff.mean / curve.sd
		p.values <- 2 * (1 - pt(abs(t.val), N.g1 - 1))
	} else {
		t.num <- (curve.mean1 - curve.mean2)
		t.den <- sqrt((N.g1 * (N.g1 - 1) * curve.sd1 ^ 2 + N.g2 * (N.g2 - 1) * curve.sd2 ^ 2) / (N.g1 + N.g2 - 2) * (1 / N.g1 + 1 / N.g2))
		t.val <- t.num / t.den
		p.values <- 2 * (1 - pt(abs(t.val), (N.g1 + N.g2 - 2)))
	}

	par(mfrow = c(1,1))
	# compute t-test at each time point
	ticks <- seq(0, max(time.all), round(max(time.all) / 10))
	plot(NULL, ,xlim = c(0, max(time.all)), ylim = c(0,1), ylab = 'Proportion of Fixations',
		xlab = 'Time', axes = FALSE, main = 'Double-Gauss Curve')
	axis(1, at = ticks)
	axis(2)
	box()
	legend('topleft', lty = 1:2, legend = groups)

	#Entries in tsmultcomp:
	#1 : Estimate of rho
	#2 : Overall Type I Error
	#3 : Total tests we will perform
	#rho.est <- ar(t.val, FALSE, order.max = 1)$ar
	rho.est <- ar1.solve(t.val)
	if(p.adj == "oleson") {
		if(paired) {
			alphastar <- find.mod.alpha(rho.est, alpha = alpha, n = N.tests, df = N.g1 - 1, cores = min(2, cores), method = "t")
		} else {
			alphastar <- find.mod.alpha(rho.est, alpha = alpha, n = N.tests, df = N.g1 + N.g2 - 2, cores = min(2, cores), method = "t")
		}
		sig <- p.values <= alphastar
	} else if(p.adj == "fdr") {
		sig <- p.adjust(p.values, "fdr") <= alpha
	} else if(p.adj == "none") {
		sig <- p.values <= alpha
	}
	
	#Make significant area yellow
	buck <- bucket(sig, time.all, ylim = c(0, .9))
	
	#Plot overall estimate of curves
	lines(time.all, curve.g1, lty = 1, lwd = 2)
	lines(time.all, curve.g2, lty = 2, lwd = 2)
	
	#Plot Confidence Interval for Group 1 curve
	lines(time.all, curve.g1 - curve.sd1 * qt(alpha / 2, N.g1 - 1), lty = 1, lwd = 1,
		col = "gray44")
	lines(time.all, curve.g1 + curve.sd1 * qt(alpha / 2, N.g1 - 1), lty = 1, lwd = 1,
		col = "gray44")
	
	#Plot Confidence Interval for Group 2 curve
	lines(time.all, curve.g2 - curve.sd2 * qt(alpha / 2, N.g2 - 1), lty = 2, lwd = 1,
		col = "gray44")
	lines(time.all, curve.g2 + curve.sd2 * qt(alpha / 2, N.g2 - 1), lty = 2, lwd = 1,
		col = "gray44")
		
	# Record confidence intervals
	curve.ci1 <- curve.ci2 <- matrix(NA, nrow = length(time.all), ncol = 4)
	curve.ci1[,1] <- curve.ci2[,1] <- time.all
	curve.ci1[,2] <- curve.g1 - curve.sd1 * qt(1 - alpha / 2, N.g1 - 1)
	curve.ci1[,3] <- curve.g1
	curve.ci1[,4] <- curve.g1 + curve.sd1 * qt(1 - alpha / 2, N.g1 - 1)
	curve.ci2[,2] <- curve.g2 - curve.sd2 * qt(1 - alpha / 2, N.g2 - 1)
	curve.ci2[,3] <- curve.g2
	curve.ci2[,4] <- curve.g2 + curve.sd2 * qt(1 - alpha / 2, N.g2 - 1)
	colnames(curve.ci1) <- colnames(curve.ci2) <- c("Time", "Lower CI", "Estimate", "Upper CI")
		
	if(!is.null(time.test)) {
		time.test <- which(time.all %in% time.test)
		
		cat("######################\n")
		cat("## Individual Tests ##\n")
		cat("######################\n")
		
		for(i in 1:length(time.test)) {
			time <- time.test[i]
			mean.1 <- curve.g1[time]
			mean.2 <- curve.g2[time]
			sd.1 <- curve.sd1[time]
			sd.2 <- curve.sd2[time]
			
			time.mean <- mean.1 - mean.2
			time.se <- sqrt((N.g1 * (N.g1 - 1) * sd.1 ^ 2 + N.g2 * (N.g2 - 1) * sd.2 ^ 2) / (N.g1 + N.g2 - 2) * (1 / N.g1 + 1 / N.g2))
			time.df <- N.g1 + N.g2 - 2
			time.t <- time.mean / time.se
			time.p <- pt(abs(time.t), time.df, lower.tail = FALSE) * 2
			pooled.sd <- sqrt((N.g1 * (N.g1 - 1) * sd.1 ^ 2 + N.g2 * (N.g2 - 1) * sd.2 ^ 2) / (N.g1 + N.g2 - 2))
			time.d <- time.mean / pooled.sd
			
			cat(paste0("Test # = ", i, " --- Time = ", time.all[time], "\n"))
			cat(paste0("Mean Diff = ", round(time.mean, 4), " --- SE = ", round(time.se, 4), "\n"))
			if(time.p < .0001) {
				cat(paste0("t = ", round(time.t, 2), " --- DF = ", round(time.df, 1), " --- p < 0.0001 \n"))
			} else {
				cat(paste0("t = ", round(time.t, 2), " --- DF = ", round(time.df, 1), " --- p = ", round(time.p, 4), "\n"))
			}
			cat(paste0("Pooled SD = ", round(pooled.sd, 4), " --- Cohen's d = ", round(time.d, 1), "\n\n"))
		}
	}
	
	if(test.params) {
		if(paired) {
			cat("######################\n")
			cat("## Parameter Tests  ##\n")
			cat("##  Paired t-test   ##\n")
			cat("######################\n")
			
			df <- N.g1 - 1
			
			mu <- mu.1 - mu.2
			mu.mean <- mean(mu)
			mu.se <- sd(mu)
			mu.t <- mu.mean / mu.se
			mu.p <- pt(abs(mu.t), df, lower.tail = FALSE) * 2
			
			ht <- ht.1 - ht.2
			ht.mean <- mean(ht)
			ht.se <- sd(ht)
			ht.t <- ht.mean / ht.se
			ht.p <- pt(abs(ht.t), df, lower.tail = FALSE) * 2
			
			s1 <- s1.1 - s1.2
			s1.mean <- mean(s1)
			s1.se <- sd(s1)
			s1.t <- s1.mean / s1.se
			s1.p <- pt(abs(s1.t), df, lower.tail = FALSE) * 2
			
			s2 <- s2.1 - s2.2
			s2.mean <- mean(s2)
			s2.se <- sd(s2)
			s2.t <- s2.mean / s2.se
			s2.p <- pt(abs(s2.t), df, lower.tail = FALSE) * 2
			
			b1 <- b1.1 - b1.2
			b1.mean <- mean(b1)
			b1.se <- sd(b1)
			b1.t <- b1.mean / b1.se
			b1.p <- pt(abs(b1.t), df, lower.tail = FALSE) * 2
			
			b2 <- b2.1 - b2.2
			b2.mean <- mean(b2)
			b2.se <- sd(b2)
			b2.t <- b2.mean / b2.se
			b2.p <- pt(abs(b2.t), df, lower.tail = FALSE) * 2
			
			cat(paste0("Mu -- Diff: ", round(mu.mean, 4), ", t: ", round(mu.t, 3), ", SE: ", round(mu.se, 3), ", df: ", df, ", p: ", round(mu.p, 4), "\n"))
			cat(paste0("Height -- Diff: ", round(ht.mean, 4), ", t: ", round(ht.t, 3), ", SE: ", round(ht.se, 3), ", df: ", df,  ", p: ", round(ht.p, 4), "\n"))
			cat(paste0("SD 1 -- Diff: ", round(s1.mean, 4), ", t: ", round(s1.t, 3), ", SE: ", round(s1.se, 3), ", df: ", df,  ", p: ", round(s1.p, 4), "\n"))
			cat(paste0("SD 2 -- Diff: ", round(s2.mean, 4), ", t: ", round(s2.t, 3), ", SE: ", round(s2.se, 3), ", df: ", df,  ", p: ", round(s2.p, 4), "\n"))
			cat(paste0("Base 1 -- Diff: ", round(b1.mean, 4), ", t: ", round(b1.t, 3), ", SE: ", round(b1.se, 3), ", df: ", df,  ", p: ", round(b1.p, 4), "\n"))
			cat(paste0("Base 2 -- Diff: ", round(b2.mean, 4), ", t: ", round(b2.t, 3), ", SE: ", round(b2.se, 3), ", df: ", df,  ", p: ", round(b2.p, 4), "\n\n"))
		} else {
			cat("######################\n")
			cat("## Parameter Tests  ##\n")
			cat("## 2 Sample t-test  ##\n")
			cat("######################\n")
			
			df <- N.g1 + N.g2 - 2
		
			mu.mean <- mean(mu.1) - mean(mu.2)
			mu.se1 <- sd(mu.1)
			mu.se2 <- sd(mu.2)
			mu.se <- sqrt((N.g1 * (N.g1 - 1) * mu.se1 ^ 2 + N.g2 * (N.g2 - 1) * mu.se2 ^ 2) / (N.g1 + N.g2 - 2) * (1 / N.g1 + 1 / N.g2))
			mu.t <- mu.mean / mu.se
			mu.p <- pt(abs(mu.t), df, lower.tail = FALSE) * 2
			
			ht.mean <- mean(ht.1) - mean(ht.2)
			ht.se1 <- sd(ht.1)
			ht.se2 <- sd(ht.2)
			ht.se <- sqrt((N.g1 * (N.g1 - 1) * ht.se1 ^ 2 + N.g2 * (N.g2 - 1) * ht.se2 ^ 2) / (N.g1 + N.g2 - 2) * (1 / N.g1 + 1 / N.g2))
			ht.t <- ht.mean / ht.se
			ht.p <- pt(abs(ht.t), df, lower.tail = FALSE) * 2
			
			s1.mean <- mean(s1.1) - mean(s1.2)
			s1.se1 <- sd(s1.1)
			s1.se2 <- sd(s1.2)
			s1.se <- sqrt((N.g1 * (N.g1 - 1) * s1.se1 ^ 2 + N.g2 * (N.g2 - 1) * s1.se2 ^ 2) / (N.g1 + N.g2 - 2) * (1 / N.g1 + 1 / N.g2))
			s1.t <- s1.mean / s1.se
			s1.p <- pt(abs(s1.t), df, lower.tail = FALSE) * 2
			
			s2.mean <- mean(s2.1) - mean(s2.2)
			s2.se1 <- sd(s2.1)
			s2.se2 <- sd(s2.2)
			s2.se <- sqrt((N.g1 * (N.g1 - 1) * s2.se1 ^ 2 + N.g2 * (N.g2 - 1) * s2.se2 ^ 2) / (N.g1 + N.g2 - 2) * (1 / N.g1 + 1 / N.g2))
			s2.t <- s2.mean / s2.se
			s2.p <- pt(abs(s2.t), df, lower.tail = FALSE) * 2
			
			b1.mean <- mean(b1.1) - mean(b1.2)
			b1.se1 <- sd(b1.1)
			b1.se2 <- sd(b1.2)
			b1.se <- sqrt((N.g1 * (N.g1 - 1) * b1.se1 ^ 2 + N.g2 * (N.g2 - 1) * b1.se2 ^ 2) / (N.g1 + N.g2 - 2) * (1 / N.g1 + 1 / N.g2))
			b1.t <- b1.mean / b1.se
			b1.p <- pt(abs(b1.t), df, lower.tail = FALSE) * 2
			
			b2.mean <- mean(b2.1) - mean(b2.2)
			b2.se1 <- sd(b2.1)
			b2.se2 <- sd(b2.2)
			b2.se <- sqrt((N.g1 * (N.g1 - 1) * b2.se1 ^ 2 + N.g2 * (N.g2 - 1) * b2.se2 ^ 2) / (N.g1 + N.g2 - 2) * (1 / N.g1 + 1 / N.g2))
			b2.t <- b2.mean / b2.se
			b2.p <- pt(abs(b2.t), df, lower.tail = FALSE) * 2
			
			cat(paste0("Mu -- Diff: ", round(mu.mean, 4), ", t: ", round(mu.t, 3), ", SE: ", round(mu.se, 3), ", df: ", df, ", p: ", round(mu.p, 4), "\n"))
			cat(paste0("Height -- Diff: ", round(ht.mean, 4), ", t: ", round(ht.t, 3), ", SE: ", round(ht.se, 3), ", df: ", df,  ", p: ", round(ht.p, 4), "\n"))
			cat(paste0("SD 1 -- Diff: ", round(s1.mean, 4), ", t: ", round(s1.t, 3), ", SE: ", round(s1.se, 3), ", df: ", df,  ", p: ", round(s1.p, 4), "\n"))
			cat(paste0("SD 2 -- Diff: ", round(s2.mean, 4), ", t: ", round(s2.t, 3), ", SE: ", round(s2.se, 3), ", df: ", df,  ", p: ", round(s2.p, 4), "\n"))
			cat(paste0("Base 1 -- Diff: ", round(b1.mean, 4), ", t: ", round(b1.t, 3), ", SE: ", round(b1.se, 3), ", df: ", df,  ", p: ", round(b1.p, 4), "\n"))
			cat(paste0("Base 2 -- Diff: ", round(b2.mean, 4), ", t: ", round(b2.t, 3), ", SE: ", round(b2.se, 3), ", df: ", df,  ", p: ", round(b2.p, 4), "\n\n"))
		}
	}
	
	params <- round(c(alpha = alpha, alpha.adj = alphastar, rho.est = rho.est), 4)
	print(list(alpha = params, significant = buck))
	invisible(list(alpha = params, significant = buck, time.all = time.all,sig = sig,
		curve.ci1 = curve.ci1, curve.ci2 = curve.ci2, curve.g1 = curve.g1, curve.g2 = curve.g2,
		curve.sd1 = curve.sd1, curve.sd2 = curve.sd2, N.g1 = N.g1, N.g2 = N.g2,
    curve1.mat = curve1.mat, curve2.mat = curve2.mat, groups = groups, seed = seed))
}