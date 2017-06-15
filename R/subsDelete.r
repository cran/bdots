subs.delete <- function(part1.list, subjs, groups) {
	if(length(subjs) != length(groups)) stop("subjs and groups should be the same length")
  
  if(any(!(groups %in% part1.list$groups))) stop("Some values in 'groups' input are not valid groups from the data set")
  groups <- ifelse(groups == part1.list$groups[1], 1, 2)
  
  if(!all((groups == 1 & subjs %in% part1.list$id.nums.g1) | (groups == 2 & subjs %in% part1.list$id.nums.g2))) stop("Some subject numbers aren't in corresponding groups")
	subs.1 <- subs.2 <- numeric(0)
  if(any(groups == 1)) subs.1 <- vapply(subjs[groups == 1], function(x) which(part1.list$id.nums.g1 == x), numeric(1))
	if(any(groups == 2)) subs.2 <- vapply(subjs[groups == 2], function(x) which(part1.list$id.nums.g2 == x), numeric(1))
	
	group.1 <- as.character(part1.list$groups[1])
	group.2 <- as.character(part1.list$groups[2])
	
	if(length(subs.1) != 0) {
		part1.list$data <- subset(part1.list$data, !(part1.list$data$Group == group.1 & part1.list$data$Subject %in% subs.1))
		part1.list$N.sub1 <- part1.list$N.sub1 - length(subs.1)
		part1.list$N.g1 <- part1.list$N.g1 - length(subs.1)
		part1.list$coef.id1 <- part1.list$coef.id1[-subs.1,]
		part1.list$coef.id3 <- part1.list$coef.id3[-subs.1,]
		part1.list$sdev.id1 <- part1.list$sdev.id1[-subs.1,]
		part1.list$sdev.id3 <- part1.list$sdev.id3[-subs.1,]
		part1.list$sigma.id1 <- part1.list$sigma.id1[-subs.1,]
		part1.list$sigma.id3 <- part1.list$sigma.id3[-subs.1,]
		part1.list$id.nums.g1 <- part1.list$id.nums.g1[-subs.1]
		part1.list$cor.1 <- part1.list$cor.1[-subs.1]
		part1.list$cor.3 <- part1.list$cor.3[-subs.1]
		part1.list$R2.g1.1 <- part1.list$R2.g1.1[-subs.1]
		part1.list$R2.g1.2 <- part1.list$R2.g1.2[-subs.1]
	}
	if(length(subs.2) != 0) {
		part1.list$data <- subset(part1.list$data, !(part1.list$data$Group == group.2 & part1.list$data$Subject %in% subs.2))
		part1.list$N.sub2 <- part1.list$N.sub2 - length(subs.2)
		part1.list$N.g2 <- part1.list$N.g2 - length(subs.2)
		part1.list$coef.id2 <- part1.list$coef.id2[-subs.2,]
		part1.list$coef.id4 <- part1.list$coef.id4[-subs.2,]
		part1.list$sdev.id2 <- part1.list$sdev.id2[-subs.2,]
		part1.list$sdev.id4 <- part1.list$sdev.id4[-subs.2,]
		part1.list$sigma.id2 <- part1.list$sigma.id2[-subs.2,]
		part1.list$sigma.id4 <- part1.list$sigma.id4[-subs.2,]
		part1.list$id.nums.g2 <- part1.list$id.nums.g2[-subs.2]
		part1.list$cor.2 <- part1.list$cor.2[-subs.2]
		part1.list$cor.4 <- part1.list$cor.4[-subs.2]
		part1.list$R2.g2.1 <- part1.list$R2.g2.1[-subs.2]
		part1.list$R2.g2.2 <- part1.list$R2.g2.2[-subs.2]
	}
	
	part1.list
}
