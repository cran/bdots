## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, include = FALSE---------------------------------------------------
library(data.table)
library(bdots)
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
            group = c("Group", "LookType"),
            curveFun = doubleGauss(concave = TRUE),
            cor = TRUE,
            numRefits = 2,
            cores = 2,
            verbose = FALSE)

refit <- brefit(fit, quickRefit = TRUE, fitCode = 5)

## -----------------------------------------------------------------------------
parDT <- coefWriteout(refit)
head(parDT)

## -----------------------------------------------------------------------------
## Subject, Group, and LookType
head(refit)

## doubleGauss pars
colnames(coef(refit))

## ----eval = FALSE-------------------------------------------------------------
# ## Save this for later using data.table::fwrite
# fwrite(parDT, file = "mypars.csv")
# parDT <- fread("mypars.csv")

## -----------------------------------------------------------------------------
new_refit <- brefit(refit, paramDT = parDT)

## -----------------------------------------------------------------------------
head(new_refit)

