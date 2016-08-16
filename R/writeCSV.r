bdots.write.csv <- function(part1.list, part2.list, file, agg = TRUE, ...) {
	if(agg) {
		group.1 <- as.character(part2.list$groups[1])
		group.2 <- as.character(part2.list$groups[2])
		
		if(part1.list$diffs) {
			agg <- with(part1.list, aggregate(data[,4], list(Time = data$Time, Curve = data$Curve, Group = data$Group), mean))

			agg.curve1 <- subset(agg, agg$Curve == 1)
			agg.curve1 <- agg.curve1[,-2]
			names(agg.curve1)[3] <- "x1"
			agg.curve2 <- subset(agg, agg$Curve == 2)
			agg.curve2 <- agg.curve2[,-2]
			names(agg.curve2)[3] <- "x2"

			agg <- merge(agg.curve1, agg.curve2)
			agg$x <- agg$x1 - agg$x2
			agg <- agg[,-c(3,4)]
			agg <- agg[order(agg$Time),]
		} else {
			agg <- with(part1.list, aggregate(data[,4], list(Time = data$Time, Group = data$Group), mean))
		}
		
		agg.1 <- subset(agg, agg$Group == part1.list$groups[1])
		agg.2 <- subset(agg, agg$Group == part1.list$groups[2])

		out <- cbind(agg.1[,3], agg.2[,3], part2.list$curve.ci1, part2.list$curve.ci2, part2.list$sig)
		out <- out[,-7]
		out <- out[,c(3, 1, 2, 4:ncol(out))]
		colnames(out)[2:10] <- c(paste(group.1, "- Observed"), paste(group.2, "- Observed"),
														paste(group.1, "- Lower CI"), paste(group.1, "- Estimate"),
														paste(group.1, "- Upper CI"), paste(group.2, "- Lower CI"),
														paste(group.2, "- Estimate"), paste(group.2, "- Upper CI"),
														"Significance")
	} else {
		
	}
	
	write.csv(out, file, ...)
}