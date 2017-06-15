subs.plot <- function(part1.list, legend.spot = "topright", ylim = NULL, groups = NA, subjs = NA, curves = NA) {
	par(mfrow = c(1,1))
	
	if(!is.na(groups) && !is.na(subjs) && length(groups) != length(subjs)) stop("Length of groups and subjs need to be equal (or left as NA to show all fits)")
  if(!is.na(groups) && any(!(groups %in% part1.list$groups))) stop("Some values in 'groups' input are not valid groups from the data set")
  if(!is.na(groups)) groups <- ifelse(groups == part1.list$groups[1], 1, 2)
  if(!is.na(groups) && !is.na(subjs) && !all((groups == 1 & subjs %in% part1.list$id.nums.g1) | (groups == 2 & subjs %in% part1.list$id.nums.g2))) stop("Some subject numbers aren't in corresponding groups")
	if(!is.na(subjs) && any(groups == 1)) subjs[groups == 1] <- vapply(subjs[groups == 1], function(x) which(part1.list$id.nums.g1 == x), numeric(1))
	if(!is.na(subjs) && any(groups == 2)) subjs[groups == 2] <- vapply(subjs[groups == 2], function(x) which(part1.list$id.nums.g2 == x), numeric(1))
  
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
	time.all   <- part1.list$time.all
	N.g1       <- part1.list$N.g1
	N.g2       <- part1.list$N.g2
	model      <- part1.list$model
	diffs      <- part1.list$diffs
	R2.g1.1    <- part1.list$R2.g1.1
	R2.g1.2    <- part1.list$R2.g1.2
	R2.g2.1    <- part1.list$R2.g2.1
	R2.g2.2    <- part1.list$R2.g2.2
	cor.1      <- part1.list$cor.1
	cor.2      <- part1.list$cor.2
	cor.3      <- part1.list$cor.3
	cor.4      <- part1.list$cor.4
	
	if(diffs && !all(is.na(groups)) && !all(is.na(curves)) && (length(curves) != length(groups))) stop("curves vector needs to be the same length as groups and subjs")
  
	if (all(is.na(subjs)) && !diffs) {
		groups <- rep(1:2, c(N.sub1, N.sub2))
		subjs <- c(seq(1, N.sub1), seq(1, N.sub2))
		curves <- rep(1, N.sub1 + N.sub2)
	} else if(all(is.na(subjs)) && diffs) {
		groups <- rep(1:2, 2 * c(N.sub1, N.sub2))
		subjs <- rep(c(seq(1, N.sub1), seq(1, N.sub2)), each = 2)
		curves <- rep(1:2, N.sub1 + N.sub2)
	} else if(diffs && !all(is.na(subjs)) && all(is.na(curves))) {
    curves <- rep(1:2, length(groups))
    groups <- rep(groups, each = 2)
    subjs <- rep(subjs, each = 2)
  }
	
	if(is.null(ylim)) {
		ylim <- c(min(data[,col]), max(data[,col]))
	}
	xticks <- seq(0, max(time.all), round(max(time.all) / 10))
	yticks <- seq(ylim[1], ylim[2], round(diff(ylim) / 10, 2))
	gray.subs <- rgb(0, 0, 0, alpha = 20, maxColorValue = 255)
	
	dgauss <- function(time, mu, ht, sig1, sig2, base1, base2) {
		(time < mu) * (exp(-1 * (time - mu) ^ 2 / (2 * sig1 ^ 2)) * (ht - base1) +
		base1) + (mu <= time) * (exp(-1 * (time - mu) ^ 2 / (2 * sig2 ^ 2)) * (ht - base2) +
		base2)
	}
	
	logist <- function(time, mini, peak, slope, cross) {
		mini + (peak-mini) / (1 + exp(4 * slope * (cross - (time)) / (peak - mini)))
	}
	
	fit.fixations <- function(model, coef) {
		if(model == "logistic") {
			mini  <- coef[1]
			peak  <- coef[2]
			slope <- coef[3]
			cross <- coef[4]
			y.fix.fitted <- logist(time.all, mini, peak, slope, cross)
		} else if (model == "dgauss") {
			mu    <- coef[1]
			ht    <- coef[2]
			sig1  <- coef[3]
			sig2  <- coef[4]
			base1 <- coef[5]
			base2 <- coef[6]
			y.fix.fitted <- dgauss(time.all, mu, ht, sig1, sig2, base1, base2)
		}
		
		y.fix.fitted
	}
	
	i <- 0
	draw <- FALSE
  while(i < length(groups)) {
    i <- i + 1
		if(groups[i] == 1) {
      subj <- subjs[i]
			id.nums <- id.nums.g1[subj]
			coef <- ifelse(rep(curves[i] == 1, ncol(coef.id1)), coef.id1[subj,], coef.id3[subj,])
			R2 <- ifelse(curves[i] == 1, R2.g1.1[subj], R2.g1.2[subj])
			cor <- ifelse(curves[i] == 1, cor.1[subj], cor.3[subj])
		} else {
      subj <- subjs[i]
			id.nums <- id.nums.g2[subj]
			coef <- ifelse(rep(curves[i] == 1, ncol(coef.id2)), coef.id2[subj,], coef.id4[subj,])
			R2 <- ifelse(curves[i] == 1, R2.g2.1[subj], R2.g2.2[subj])
			cor <- ifelse(curves[i] == 1, cor.2[subj], cor.4[subj])
		}
		
		if(diffs) {
			y1id <- subset(data, data$Subject == id.nums & data$Group == part1.list$groups[groups[i]] & data$Curve == curves[i])
		} else {
			y1id <- subset(data, data$Subject == id.nums & data$Group == part1.list$groups[groups[i]])
		}
		y.fix <- y1id[,col]
		
		y.fix.fitted <- fit.fixations(model, coef)
		
		plot(time.all, y.fix, type = "l", lwd = 2, xlab = "Time", ylab = "Proportion of Fixations",
				ylim = ylim, axes = FALSE)
		axis(1, at = xticks)
		axis(2, at = yticks)
		box()
		abline(v = xticks, h = yticks, col = gray.subs)
		if(diffs) {
			title(paste0("Group = ", part1.list$groups[groups[i]], ", Subject = ", id.nums, ", Curve = ", curves[i], ", R2 = ", round(R2, 3), ", Cor = ", ifelse(cor == 1, "AR", "Non-AR")))
		} else {
			title(paste0("Group = ", part1.list$groups[groups[i]], ", Subject = ", id.nums, ", Curve = 1", ", R2 = ", round(R2, 3), ", Cor = ", ifelse(cor == 1, "AR", "Non-AR")))
		}
		if(model == "logistic") {
			title(sub = paste0("Mini = ", round(coef[1], 2), ", Peak = ", round(coef[2], 2), ", Slope = ", round(coef[3], 4), ", Cross = ", round(coef[4])))
		} else if(model == "dgauss") {
			title(sub = paste0("Mu = ", round(coef[1]), ", Height = ", round(coef[2], 2), ", Sigma1 = ", round(coef[3]), ", Sigma2 = ", round(coef[4]),
													", Base1 = ", round(coef[5], 3), ", Base2 = ", round(coef[6], 3)))
		}
		lines(time.all, y.fix.fitted, lty = 2, lwd = 2, col = "red")
		if(draw) {
			y.fix.draw <- fit.fixations(model, draw.coef)
			lines(time.all, y.fix.draw, lty = 3, lwd = 2, col = "blue")
			draw <- FALSE
		}
		legend(legend.spot, lty = 1:2, col = c("black", "red"), legend = c("Observed", "Fitted"))
		
		
		
		draw <- FALSE
		read.good <- FALSE
		while(!read.good) {
			read <- readline("Press 'Enter' to go to next plot, type 'b' to go back 1 plot")
			
			read.good <- TRUE
			if(read == "back" || read == "b") {
				i <- ifelse(i > 1, i - 2, 0)
			} else if(grepl("^d(raw)? ", read)) {
				draw.spots <- gregexpr(" ", read)
				draw.len <- as.numeric(attr(draw.spots[[1]], "match.length"))
				draw.spots <- as.numeric(draw.spots[[1]])
				
				if(length(draw.spots) != 4 && model == "logistic") {
					warning("Model is the 4-parameter logistic,\nplease specify 4 parameters with the following form:\n\"d B P S X\"\n")
					read.good <- FALSE
				} else if(length(draw.spots) !=6 && model == "dgauss") {
					warning("Model is the double gauss,\nplease specify 6 parameters with the following form:\n\"B1 B2 S1 S2 Mu P\"")
					read.good <- FALSE
				} else {
					draw <- TRUE
					i <- i - 1
					draw.coef <- apply(cbind(draw.spots, draw.len), 1, function(x) substr(read, x[1] + 1, x[1] + x[2]))
				}
			}
		}
	}
}