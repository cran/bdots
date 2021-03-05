## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(bdots)

## -----------------------------------------------------------------------------
fit <- bdotsFit(data = cohort_unrelated,
                subject = "Subject",
                time = "Time",
                y = "Fixations",
                group = c("DB_cond", "LookType"),
                curveType = doubleGauss(concave = TRUE),
                cores = 2)

## -----------------------------------------------------------------------------
doubleGauss

## -----------------------------------------------------------------------------
## Return a unique subject/group permutation
dat <- cohort_unrelated[Subject == 1 & DB_cond == 50 & LookType == "Cohort", ]
dat

## -----------------------------------------------------------------------------
## See return value
doubleGauss(dat = dat, y = "Fixations", time = "Time", concave = TRUE)

## -----------------------------------------------------------------------------
doubleGauss2 <- function (dat, params = NULL, concave = TRUE, ...) {
    dgaussPars <- function(dat, conc = concave) {
        time <- dat[["Time"]]
        y <- dat[["Fixations"]]
        mu <- ifelse(conc, time[which.max(y)], time[which.min(y)])
        ht <- ifelse(conc, max(y), min(y))
        base1 <- ifelse(conc, min(y[time < mu]), max(y[time < 
            mu]))
        base2 <- ifelse(conc, min(y[time > mu]), max(y[time > 
            mu]))
        y1 <- y - base1
        y1 <- rev(y1[time <= mu])
        time1 <- rev(time[time <= mu])
        totalY1 <- sum(y1)
        sigma1 <- mu - time1[which.min(abs((pnorm(1) - pnorm(-1)) * 
            totalY1 - cumsum(y1)))]
        y2 <- y - base2
        y2 <- y2[time >= mu]
        time2 <- time[time >= mu]
        totalY2 <- sum(y2)
        sigma2 <- time2[which.min(abs((pnorm(1) - pnorm(-1)) * 
            totalY2 - cumsum(y2)))] - mu
        return(c(mu = mu, ht = ht, sig1 = sigma1, sig2 = sigma2, 
            base1 = base1, base2 = base2))
    }
    if (is.null(params)) {
        params <- dgaussPars(dat, concave)
    }
    else {
        if (length(params) != 6) 
            stop("doubleGauss requires 6 parameters be specified for refitting")
        if (!all(names(params) %in% c("mu", "ht", "sig1", "sig2", 
            "base1", "base2"))) {
            stop("doubleGauss parameters for refitting must be correctly labeled")
        }
    }

    ## Here, we use Fixations and Time directly
    ff <- bquote(Fixations ~ (Time < mu) * (exp(-1 * (Time - 
        mu)^2/(2 * sig1^2)) * (ht - base1) + base1) + (mu <= 
        Time) * (exp(-1 * (Time - mu)^2/(2 * sig2^2)) * 
        (ht - base2) + base2))
    return(list(formula = ff, params = params))
}

same_fit_different_day <- bdotsFit(data = cohort_unrelated,
                                   subject = "Subject",
                                   time = "Time",
                                   y = "Fixations",
                                   group = c("DB_cond", "LookType"),
                                   curveType = doubleGauss2(concave = TRUE),
                                   cores = 2)

## -----------------------------------------------------------------------------
## Original fit
head(coef(fit))

## "New" fit
head(coef(same_fit_different_day))

