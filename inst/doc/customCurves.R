## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, include = FALSE---------------------------------------------------
library(bdots)
library(data.table)
# Make smaller for cran
cohort_unrelated$Subject <- as.numeric(cohort_unrelated$Subject)
cohort_unrelated <- as.data.table(cohort_unrelated)
cohort_unrelated <- cohort_unrelated[Subject < 10, ]

## -----------------------------------------------------------------------------
library(bdots)

fit <- bfit(data = cohort_unrelated,
            subject = "Subject",
            time = "Time",
            y = "Fixations",
            group = c("DB_cond", "LookType"),
            curveFun = doubleGauss(concave = TRUE),
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
doubleGauss2 <- function (dat, y, time, params = NULL, concave = TRUE, ...) {
  
  if (is.null(params)) {
    ## Instead of defining our own, just reuse the one in bdots
    params <- bdots:::dgaussPars(dat, y, time, concave)
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
    ff <- bquote(Fixations ~ (Time < mu) * (exp(-1 * (Time - mu)^2 / 
                  (2 * sig1^2)) * (ht - base1) + base1) + (mu <= Time) * 
                  (exp(-1 * (Time - mu)^2/(2 * sig2^2)) * (ht - base2) + base2))
    return(list(formula = ff, params = params))
}

same_fit_different_day <- bfit(data = cohort_unrelated,
                               subject = "Subject",
                               time = "Time",
                               y = "Fixations",
                               group = c("DB_cond", "LookType"),
                               curveFun = doubleGauss2(concave = TRUE),
                               cores = 2)

## -----------------------------------------------------------------------------
## Original fit
head(coef(fit))

## "New" fit
head(coef(same_fit_different_day))

