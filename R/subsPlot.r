subs.plot <- function(part1.list, legend.spot = "topright") {
	par(mfrow = c(1,1))
	
	data       <- part1.list$data
  col        <- part1.list$col
  N.time     <- part1.list$N.time
  N.sub1     <- part1.list$N.sub1
  N.sub2     <- part1.list$N.sub2
  coef.id1   <- part1.list$coef.id1
  coef.id2   <- part1.list$coef.id2
	coef.id3   <- part1.list$coef.id3
  coef.id4   <- part1.list$coef.id4
  id.nums.g1 <- part1.list$id.nums.g1
  id.nums.g2 <- part1.list$id.nums.g2
  groups     <- part1.list$groups
  time.all   <- part1.list$time.all
  N.g1       <- part1.list$N.g1
  N.g2       <- part1.list$N.g2
	model      <- part1.list$model
	diffs      <- part1.list$diffs
	
	dgauss <- function(time, mu, ht, sig1, sig2, base1, base2) {
		(time < mu) * (exp(-1 * (time - mu) ^ 2 / (2 * sig1 ^ 2)) * (ht - base1) +
		base1) + (mu <= time) * (exp(-1 * (time - mu) ^ 2 / (2 * sig2 ^ 2)) * (ht - base2) +
		base2)
	}
	
	logist <- function(time, mini, peak, slope, cross) {
		mini + (peak-mini) / (1 + exp(4 * slope * (cross - (time)) / (peak - mini)))
	}
	
	for(i in 1:N.g1) {
		if(diffs) {
			y1id <- subset(data, data$Subject == id.nums.g1[i] & data$Group == groups[1] & data$Curve == 1)
		} else {
			y1id <- subset(data, data$Subject == id.nums.g1[i] & data$Group == groups[1])
		}
		y.fix <- y1id[,col]
		
		if(model == "logistic") {
			mini  <- coef.id1[i,1]
			peak  <- coef.id1[i,2]
			slope <- coef.id1[i,3]
			cross <- coef.id1[i,4]
			y.fix.fitted <- logist(time.all, mini, peak, slope, cross)
		} else if (model == "dgauss") {
			mu    <- coef.id1[i,1]
			ht    <- coef.id1[i,2]
			sig1  <- coef.id1[i,3]
			sig2  <- coef.id1[i,4]
			base1 <- coef.id1[i,5]
			base2 <- coef.id1[i,6]
			y.fix.fitted <- dgauss(time.all, mu, ht, sig1, sig2, base1, base2)
		}

		plot(time.all, y.fix, type = "l", lwd = 2, xlab = "Time", ylab = "Proportion of Fixations")
		if(diffs) {
			title(paste0("Group = ", groups[1], ", Subject = ", id.nums.g1[i], ", Dataset # =", i, ", Curve = 1"))
		} else {
			title(paste0("Group = ", groups[1], ", Subject = ", id.nums.g1[i], ", Dataset # =", i))
		}
		lines(time.all, y.fix.fitted, lty = 2, lwd = 2, col = "red")
		legend(legend.spot, lty = 1:2, col = c("black", "red"),
			legend = c("Observed", "Fitted"))
		
		readline("Press 'Enter' to go to next plot")
		
		if(diffs) {
			y1id <- subset(data, data$Subject == id.nums.g1[i] & data$Group == groups[1] & data$Curve == 2)
			y.fix <- y1id[,col]
			
			if(model == "logistic") {
				mini  <- coef.id3[i,1]
				peak  <- coef.id3[i,2]
				slope <- coef.id3[i,3]
				cross <- coef.id3[i,4]
				y.fix.fitted <- logist(time.all, mini, peak, slope, cross)
			} else if (model == "dgauss") {
				mu    <- coef.id3[i,1]
				ht    <- coef.id3[i,2]
				sig1  <- coef.id3[i,3]
				sig2  <- coef.id3[i,4]
				base1 <- coef.id3[i,5]
				base2 <- coef.id3[i,6]
				y.fix.fitted <- dgauss(time.all, mu, ht, sig1, sig2, base1, base2)
			}

			plot(time.all, y.fix, type = "l", lwd = 2, xlab = "Time", ylab = "Proportion of Fixations")
			title(paste0("Group = ", groups[1], ", Subject = ", id.nums.g1[i], ", Dataset # =", i, ", Curve = 2"))
			lines(time.all, y.fix.fitted, lty = 2, lwd = 2, col = "red")
			legend(legend.spot, lty = 1:2, col = c("black", "red"),
				legend = c("Observed", "Fitted"))
				
			readline("Press 'Enter' to go to next plot")
		}
	}
	
	for(i in 1:N.g2) {
		if(diffs) {
			y1id <- subset(data, data$Subject == id.nums.g2[i] & data$Group == groups[2] & data$Curve == 1)
		} else {
			y1id <- subset(data, data$Subject == id.nums.g2[i] & data$Group == groups[2])
		}
		
		y.fix <- y1id[,col]
		
		if(model == "logistic") {
			mini  <- coef.id2[i,1]
			peak  <- coef.id2[i,2]
			slope <- coef.id2[i,3]
			cross <- coef.id2[i,4]
			y.fix.fitted <- logist(time.all, mini, peak, slope, cross)
		} else if (model == "dgauss") {
			mu    <- coef.id2[i,1]
			ht    <- coef.id2[i,2]
			sig1  <- coef.id2[i,3]
			sig2  <- coef.id2[i,4]
			base1 <- coef.id2[i,5]
			base2 <- coef.id2[i,6]
			y.fix.fitted <- dgauss(time.all, mu, ht, sig1, sig2, base1, base2)
		}
		
		plot(time.all, y.fix, type = "l", lwd = 2, xlab = "Time", ylab = "Proportion of Fixations")
		if(diffs) {
			title(paste0("Group = ", groups[2], ", Subject = ", id.nums.g2[i], ", Dataset # =", i, ", Curve = 1"))
		} else {
			title(paste0("Group = ", groups[2], ", Subject = ", id.nums.g2[i], ", Dataset # =", i))
		}
		lines(time.all, y.fix.fitted, lty = 2, lwd = 2, col = "red")
		legend(legend.spot, lty = 1:2, col = c("black", "red"),
			legend = c("Observed", "Fitted"))
		
		readline("Press 'Enter' to go to next plot")
		
		if(diffs) {
			y1id <- subset(data, data$Subject == id.nums.g2[i] & data$Group == groups[2] & data$Curve == 2)
			
			y.fix <- y1id[,col]
			
			if(model == "logistic") {
				mini  <- coef.id4[i,1]
				peak  <- coef.id4[i,2]
				slope <- coef.id4[i,3]
				cross <- coef.id4[i,4]
				y.fix.fitted <- logist(time.all, mini, peak, slope, cross)
			} else if (model == "dgauss") {
				mu    <- coef.id4[i,1]
				ht    <- coef.id4[i,2]
				sig1  <- coef.id4[i,3]
				sig2  <- coef.id4[i,4]
				base1 <- coef.id4[i,5]
				base2 <- coef.id4[i,6]
				y.fix.fitted <- dgauss(time.all, mu, ht, sig1, sig2, base1, base2)
			}
			
			plot(time.all, y.fix, type = "l", lwd = 2, xlab = "Time", ylab = "Proportion of Fixations")
			title(paste0("Group = ", groups[2], ", Subject = ", id.nums.g2[i], ", Dataset # =", i, ", Curve = 2"))
			lines(time.all, y.fix.fitted, lty = 2, lwd = 2, col = "red")
			legend(legend.spot, lty = 1:2, col = c("black", "red"),
				legend = c("Observed", "Fitted"))
			
			readline("Press 'Enter' to go to next plot")
		}
	}
}
