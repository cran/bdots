## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=FALSE-----------------------------------------------------
library(data.table)
library(bdots)

## ----include=FALSE------------------------------------------------------------
# Make smaller for cran
cohort_unrelated$Subject <- as.numeric(cohort_unrelated$Subject)
cohort_unrelated <- as.data.table(cohort_unrelated)
cohort_unrelated <- cohort_unrelated[Subject < 10, ]

## -----------------------------------------------------------------------------
head(cohort_unrelated)

## -----------------------------------------------------------------------------
fit <- bfit(data = cohort_unrelated,
            subject = "Subject",
            time = "Time",
            y = "Fixations",
            group = c("DB_cond", "LookType"),
            curveFun = doubleGauss(concave = TRUE),
            cores = 2)

## -----------------------------------------------------------------------------
head(coef(fit))

head(coef(fit[DB_cond == 50, ]))

## ----fig.align='center', fig.width = 8, fig.height=6--------------------------
plot(fit[1:4, ])

## ----eval = FALSE-------------------------------------------------------------
# ## Quickly auto-refit (not run)
# refit <- brefit(fit, fitCode = 1L, quickRefit = TRUE)
# 
# ## Manual refit (not run)
# refit <- brefit(fit, fitCode = 1L)

## -----------------------------------------------------------------------------
table(fit$fitCode)

## Remove all failed curve fits
refit <- bdRemove(fit, fitCode = 6L)

table(refit$fitCode)

## ----eval = FALSE-------------------------------------------------------------
# ## Only one grouping variable in dataset, take bootstrapped difference
# Outcome ~ Group1(value1, value2)
# 
# ## More than one grouping variable in difference, must specify unique value
# Outcome ~ Group1(value1, value2) + Group2(value3)

## ----eval = FALSE-------------------------------------------------------------
# ## Must add LookType(Cohort) to specify
# Fixations ~ DB_cond(50, 65) + LookType(Cohort)

## ----eval = FALSE-------------------------------------------------------------
# ## Difference of difference. Here, outer difference is Group1, inner is Group2
# diffs(Outcome, Group2(value3, value4)) ~ Group1(value1, value2)
# 
# ## Same as above if three or more grouping variables
# diffs(Outcome, Group2(value3, value4)) ~ Group1(value1, value2) + Group3(value5)

## ----eval = FALSE-------------------------------------------------------------
# diffs(Fixations, DB_cond(50, 65)) ~ LookType(Cohort, Unrelated_Cohort)

## -----------------------------------------------------------------------------
boot1 <- bboot(formula = Fixation ~ DB_cond(50, 65) + LookType(Cohort),
               bdObj = refit,
               Niter = 1000,
               alpha = 0.05,
               padj = "oleson",
               cores = 2)

boot2 <- bboot(formula = diffs(Fixation, LookType(Cohort, Unrelated_Cohort)) ~ DB_cond(50, 65),
               bdObj = refit,
               Niter = 1000,
               alpha = 0.05,
               padj = "oleson",
               cores = 2)

## ----fig.align='center', fig.width = 8, fig.height=6--------------------------
summary(boot1)

plot(boot1)

