logistic.refit <- function(part1.list, subj = NULL, group = NULL, curves = NULL,
	params = NULL, cor = NULL, rho.0 = NULL, info.matrix = NULL, get.cov.only = FALSE) {
	
	#info.matrix:
	#	1st column: subject #
	# 2nd column: groups #
	# 3rd column: curve #
	# 4th column: params -- mini
	# 5th column: params -- peak
	# 6th column: params -- slope
	# 7th column: params -- cross
	# 8th column: rho.0, if column is available
	if(!is.null(info.matrix)) {
		if(any(is.na(info.matrix))) stop("Can't have any NA or NaN values in info.matrix")
		if(!is.matrix(info.matrix)) stop("info.matrix should be a matrix")
		if(ncol(info.matrix) != 7)
			stop("info.matrix should have 7 columns:
			subject, groups, curves, mini, peak, slope, cross
			If diffs=FALSE, curves column is not used, so just set to any number")
		subj <- as.numeric(info.matrix[,1])
		group <- info.matrix[,2]
		curves <- as.numeric(info.matrix[,3])
		params <- list()
		for(i in 1:nrow(info.matrix)) {
			params[[i]] <- as.numeric(info.matrix[i, 4:7])
		}
		if(ncol(info.matrix) == 8) rho.0 <- as.numeric(info.matrix[,8])
	}
	
	if(length(subj) == 0 || length(group) == 0) stop("Need entry for 'subj' and 'group'")
	if(length(subj) != length(group)) stop("Length of 'subj' and 'group' need to be the same")
	if(!is.null(params) && !is.list(params)) stop("'params' needs to be a list")
	if(!is.null(params) && is.list(params) && length(params) != length(subj))
		stop("Length of 'params' should be equal to length of 'subj' and 'group'")
	if(!is.null(cor) && !is.logical(cor)) stop("'cor' should be a vector of logicals")
	if(!is.null(cor) && length(cor) != length(subj))
		stop("Length of 'cor' should be equal to length of 'subj' and 'group'")
	if(is.null(cor)) cor <- rep(part1.list$cor, length(subj))
  
  if(any(!(group %in% part1.list$groups))) stop("Some values in 'group' input are not valid groups from the data set")
  group <- ifelse(group == part1.list$groups[1], 1, 2)
  
  if(!all((group == 1 & subj %in% part1.list$id.nums.g1) | (group == 2 & subj %in% part1.list$id.nums.g2))) stop("Some subject numbers aren't in corresponding group")
  if(any(group == 1)) subj[group == 1] <- vapply(subj[group == 1], function(x) which(part1.list$id.nums.g1 == x), numeric(1))
	if(any(group == 2)) subj[group == 2] <- vapply(subj[group == 2], function(x) which(part1.list$id.nums.g2 == x), numeric(1))
	subj <- as.numeric(subj)
	
	diffs <- part1.list$diffs
  if(diffs && length(curves) != length(subj)) stop("Length of 'curves' and 'subj' need to be the same")
	
	data <- part1.list$data
	id.nums.g1 <- part1.list$id.nums.g1
	id.nums.g2 <- part1.list$id.nums.g2
	time.all <- part1.list$time.all
	y.fix <- part1.list$y.fix
	if(is.null(rho.0)) rho.0 <- part1.list$rho.0
	if(length(rho.0) == 1) rho.0 <- rep(rho.0, length(subj))
	
	coef.id1 <- part1.list$coef.id1
	sdev.id1 <- part1.list$sdev.id1
	sdev.id1.cov <- part1.list$sdev.id1.cov
	sigma.id1 <- part1.list$sigma.id1
	
	coef.id2 <- part1.list$coef.id2
	sdev.id2 <- part1.list$sdev.id2
	sdev.id2.cov <- part1.list$sdev.id2.cov
	sigma.id2 <- part1.list$sigma.id2
	
	coef.id3 <- part1.list$coef.id3
	sdev.id3 <- part1.list$sdev.id3
	sdev.id3.cov <- part1.list$sdev.id3.cov
	sigma.id3 <- part1.list$sigma.id3
	
	coef.id4 <- part1.list$coef.id4
	sdev.id4 <- part1.list$sdev.id4
	sdev.id4.cov <- part1.list$sdev.id4.cov
	sigma.id4 <- part1.list$sigma.id4
	
	R2.g1.1 <- part1.list$R2.g1.1
	R2.g1.2 <- part1.list$R2.g1.2
	R2.g2.1 <- part1.list$R2.g2.1
	R2.g2.2 <- part1.list$R2.g2.2
	
	cor.1 <- part1.list$cor.1
	cor.2 <- part1.list$cor.2
	cor.3 <- part1.list$cor.3
	cor.4 <- part1.list$cor.4
	
  col <- part1.list$col
  groups <- part1.list$groups
	
	if(all(group %in% groups)) {
		group <- apply(matrix(group), 1, function(x) which(x == groups))
	}
	
	curve.f <- function(mini,peak,slope,cross,t)
		mini + (peak - mini) / (1 + exp(4 * slope * (cross - t) / (peak-mini)))
		
	fail <- 0
		
	for(i in 1:length(group)) {
		if(group[i] == 1) {
			if(diffs) {
				y1id <- subset(data, data$Subject == id.nums.g1[subj[i]] & data$Group == groups[1] & data$Curve == curves[i])
			} else {
				y1id <- subset(data, data$Subject == id.nums.g1[subj[i]] & data$Group == groups[1])
			}
			
			y.fix <- y1id[,col]

			fit.curve <- est.logistic.curve(time.all, y.fix, rho.0, params = params[[i]], cor = cor[i], get.cov.only = get.cov.only)
				
			if(diffs && curves[i] == 2) {
				old.R2 <- R2.g1.2[subj[i]]
				old.cor <- cor.3[subj[i]]
				if(is.null(fit.curve$fit)) {
					coef.id3[subj[i],] <- rep(NA, 4)
					sdev.id3.cov[[subj[i]]] <- matrix(NA, nrow = 4, ncol = 4)
					sigma.id3[subj[i],] <- NA
					cor.3[subj[i]] <- NA
					R2.g1.2[subj[i]] <- NA
					fail <- fail + 1
				} else {
					cor.3[subj[i]] <- fit.curve$cor
					fit.curve <- fit.curve$fit
					coef.id3[subj[i],] <- coef(fit.curve)
					sdev.id3.cov[[subj[i]]] <- fit.curve$varBeta
					sigma.id3[subj[i],] <- fit.curve$sigma
					
					SSY <- sum((y.fix - mean(y.fix)) ^ 2)
					y.fit <- curve.f(coef(fit.curve)[1], coef(fit.curve)[2], coef(fit.curve)[3],
						coef(fit.curve)[4], time.all)
					y.err <- y.fit - y.fix
					SSE <- sum(y.err ^ 2)
					R2.g1.2[subj[i]] <- 1 - SSE / SSY
				}
				cat("Subject = ", as.character(id.nums.g1[subj[i]]),
						", Group = ", as.character(groups[group[i]]),
						", Curve = ", curves[i],
						", Old R2 = ", round(old.R2, 3),
						", New R2 = ", round(R2.g1.2[subj[i]], 3),
						"\n\t Old AR1 = ", as.logical(old.cor),
						", New AR1 = ", as.logical(cor.3[subj[i]]),
						"\n")
			} else {
				old.R2 <- R2.g1.1[subj[i]]
				old.cor <- cor.1[subj[i]]
				if(is.null(fit.curve$fit)) {
					coef.id1[subj[i],] <- rep(NA, 4)
					sdev.id1.cov[[subj[i]]] <- matrix(NA, nrow = 4, ncol = 4)
					sigma.id1[subj[i],] <- NA
					cor.1[subj[i]] <- NA
					R2.g1.1[subj[i]] <- NA
					fail <- fail + 1
				} else {
					cor.1[subj[i]] <- fit.curve$cor
					fit.curve <- fit.curve$fit
					coef.id1[subj[i],] <- coef(fit.curve)
					sdev.id1.cov[[subj[i]]] <- fit.curve$varBeta
					sigma.id1[subj[i],] <- fit.curve$sigma
					
					SSY <- sum((y.fix - mean(y.fix)) ^ 2)
					y.fit <- curve.f(coef(fit.curve)[1], coef(fit.curve)[2], coef(fit.curve)[3],
						coef(fit.curve)[4], time.all)
					y.err <- y.fit - y.fix
					SSE <- sum(y.err ^ 2)
					R2.g1.1[subj[i]] <- 1 - SSE / SSY
				}
				cat("Subject = ", as.character(id.nums.g1[subj[i]]),
						", Group = ", as.character(groups[group[i]]),
						", Curve = ", curves[i],
						", Old R2 = ", round(old.R2, 3),
						", New R2 = ", round(R2.g1.1[subj[i]], 3),
						"\n\t Old AR1 = ", as.logical(old.cor),
						", New AR1 = ", as.logical(cor.1[subj[i]]),
						"\n")
		}} else {
			if(diffs) {
				y1id <- subset(data, data$Subject == id.nums.g2[subj[i]] & data$Group == groups[2] & data$Curve == curves[i])
			} else {
				y1id <- subset(data, data$Subject == id.nums.g2[subj[i]] & data$Group == groups[2])
			}

			y.fix <- y1id[,col]
			fit.curve <- est.logistic.curve(time.all, y.fix, rho.0, params = params[[i]], cor = cor[i], get.cov.only = get.cov.only)
				
			if(diffs && curves[i] == 2) {
				old.R2 <- R2.g2.2[subj[i]]
				old.cor <- cor.4[subj[i]]
				if(is.null(fit.curve$fit)) {
					coef.id4[subj[i],] <- rep(NA, 4)
					sdev.id4.cov[[subj[i]]] <- matrix(NA, nrow = 4, ncol = 4)
					sigma.id4[subj[i],] <- NA
					cor.4[subj[i]] <- NA
					R2.g2.2[subj[i]] <- NA
					fail <- fail + 1
				} else {
					cor.4[subj[i]] <- fit.curve$cor
					fit.curve <- fit.curve$fit
					coef.id4[subj[i],] <- coef(fit.curve)
					sdev.id4.cov[[subj[i]]] <- fit.curve$varBeta
					sigma.id4[subj[i],] <- fit.curve$sigma
					
					SSY <- sum((y.fix - mean(y.fix)) ^ 2)
					y.fit <- curve.f(coef(fit.curve)[1], coef(fit.curve)[2], coef(fit.curve)[3],
						coef(fit.curve)[4], time.all)
					y.err <- y.fit - y.fix
					SSE <- sum(y.err ^ 2)
					R2.g2.2[subj[i]] <- 1 - SSE / SSY
				}
				cat("Subject = ", as.character(id.nums.g2[subj[i]]),
						", Group = ", as.character(groups[group[i]]),
						", Curve = ", curves[i],
						", Old R2 = ", round(old.R2, 3),
						", New R2 = ", round(R2.g2.2[subj[i]], 3),
						"\n\t Old AR1 = ", as.logical(old.cor),
						", New AR1 = ", as.logical(cor.4[subj[i]]),
						"\n")
			} else {
				old.R2 <- R2.g2.1[subj[i]]
				old.cor <- cor.2[subj[i]]
				if(is.null(fit.curve$fit)) {
					coef.id2[subj[i],] <- rep(NA, 4)
					sdev.id2.cov[[subj[i]]] <- matrix(NA, nrow = 4, ncol = 4)
					sigma.id2[subj[i],] <- NA
					cor.2[subj[i]] <- NA
					R2.g2.1[subj[i]] <- NA
					fail <- fail + 1
				} else {
					cor.2[subj[i]] <- fit.curve$cor
					fit.curve <- fit.curve$fit
					coef.id2[subj[i],] <- coef(fit.curve)
					sdev.id2.cov[[subj[i]]] <- fit.curve$varBeta
					sigma.id2[subj[i],] <- fit.curve$sigma
					
					SSY <- sum((y.fix - mean(y.fix)) ^ 2)
					y.fit <- curve.f(coef(fit.curve)[1], coef(fit.curve)[2], coef(fit.curve)[3],
						coef(fit.curve)[4], time.all)
					y.err <- y.fit - y.fix
					SSE <- sum(y.err ^ 2)
					R2.g2.1[subj[i]] <- 1 - SSE / SSY
				}
				cat("Subject = ", as.character(id.nums.g2[subj[i]]),
						", Group = ", as.character(groups[group[i]]),
						", Curve = ", curves[i],
						", Old R2 = ", round(old.R2, 3),
						", New R2 = ", round(R2.g2.1[subj[i]], 3),
						"\n\t Old AR1 = ", as.logical(old.cor),
						", New AR1 = ", as.logical(cor.2[subj[i]]),
						"\n")
      }
		}
	}
	
	part1.list$coef.id1 <- coef.id1 
	part1.list$sdev.id1.cov <- sdev.id1.cov
	part1.list$sigma.id1 <- sigma.id1 
	part1.list$coef.id2 <- coef.id2 
	part1.list$sdev.id2.cov <- sdev.id2.cov
	part1.list$sigma.id2 <- sigma.id2
	part1.list$coef.id3 <- coef.id3 
	part1.list$sdev.id3.cov <- sdev.id3.cov
	part1.list$sigma.id3 <- sigma.id3
	part1.list$coef.id4 <- coef.id4 
	part1.list$sdev.id4.cov <- sdev.id4.cov
	part1.list$sigma.id4 <- sigma.id4
	
	part1.list$R2.g1.1 <- R2.g1.1 
	part1.list$R2.g1.2 <- R2.g1.2 
	part1.list$R2.g2.1 <- R2.g2.1 
	part1.list$R2.g2.2 <- R2.g2.2 
	
	part1.list$cor.1 <- cor.1
	part1.list$cor.2 <- cor.2
	part1.list$cor.3 <- cor.3
	part1.list$cor.4 <- cor.4
	
	if(fail > 0) warning(paste("Unable to find fits for", fail, "curves"))
	
	part1.list
}