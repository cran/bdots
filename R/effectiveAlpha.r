effectiveAlpha <- function(rho, k, n, df, method = "mvt") {
	if(method == "mvt") {
		out <- effectiveAlpha.mvT(rho, k, n, df)
	} else if(method == "mvn") {
		out <- effectiveAlpha.mvNorm(rho, k, n)
	} else {
		out <- effectiveAlpha.normApprox(rho, k, n)
	}
	
	out
}