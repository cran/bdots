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
                group = c("Group", "LookType"),
                curveType = doubleGauss(concave = TRUE),
                cor = TRUE,
                numRefits = 2,
                cores = 2,
                verbose = FALSE)

refit <- bdotsRefit(fit, quickRefit = TRUE, fitCode = 5)

## -----------------------------------------------------------------------------
parDT <- coefWriteout(refit)
head(parDT)

## -----------------------------------------------------------------------------
## Subject, Group, and LookType
head(refit)

## doubleGauss pars
colnames(coef(refit))

## ---- eval = FALSE------------------------------------------------------------
#  ## Save this for later using data.table::fwrite
#  fwrite(parDT, file = "mypars.csv")
#  parDT <- fread("mypars.csv")

## -----------------------------------------------------------------------------
new_refit <- bdotsRefit(refit, paramDT = parDT)

## -----------------------------------------------------------------------------
head(new_refit)

