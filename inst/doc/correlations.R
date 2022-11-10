## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- include=FALSE-----------------------------------------------------------
library(bdots)
# Make smaller for cran
cohort_unrelated$Subject <- as.numeric(cohort_unrelated$Subject)
cohort_unrelated <- as.data.table(cohort_unrelated)
cohort_unrelated <- cohort_unrelated[Subject < 10, ]

## -----------------------------------------------------------------------------
library(bdots)
library(data.table)

## Let's work with cohort_unrelated dataset, as it has multiple groups
dat <- as.data.table(cohort_unrelated)

## And add a fixed value for which we want to find a correlation
dat[, val := rnorm(1), by = Subject]

head(dat)

## -----------------------------------------------------------------------------
## Create regular fit in bdots
fit <- bdotsFit(data = dat,
                subject = "Subject",
                time = "Time",
                group = c("LookType", "Group"),
                y = "Fixations", curveType = doubleGauss2(),
                cores = 2)

## -----------------------------------------------------------------------------
## Returns a data.table of class bdotsCorrObj
corr_ci <- bdotsCorr(fit, val = "val", ciBands = TRUE)
head(corr_ci)

## Same, without confidence intervals
corr_noci <- bdotsCorr(fit, val = "val")
head(corr_noci)

## ---- fig.align='center', fig.width = 6, fig.height=6-------------------------
## Default is no bands
plot(corr_ci)

## Try again with bands
plot(corr_ci, ciBands = TRUE)

## Narrow in on a particular window
plot(corr_ci, window = c(750, 1500))

## ---- fig.align='center', fig.width = 6, fig.height=4-------------------------
plot(corr_ci[Group2 == "50", ])

