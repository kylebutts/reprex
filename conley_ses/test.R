setwd("~/Documents/reprex/conley_ses/")
library(Rcpp)
library(RcppArmadillo)
library(tidyverse)
library(lfe)

source("helper-conley.R")

# Test XeeX_spatial ------------------------------------------------------------
# From ARE 212 Section 10: Non-Standard Standard Errors II
# and https://github.com/edrubin/ARE212

data(quakes) 
m <- felm(depth ~ 1 + mag | 0 | 0 | lat+long, data = quakes, keepCX = TRUE)

M <- cbind(quakes$lat, quakes$long)
X <- as.matrix(cbind(m$cX))
e <- resid(m)
dist_cutoff <- 100

meat <- XeeX_spatial(M=M, cutoff= 100, X=X, e=e, k=2, kernel = "uniform")

sqrt(diag(solve(t(X) %*% X) %*% meat %*% solve(t(X) %*% X)))

## Should be ~ this (they use worse distance formula):
#  (Intercept)  mag
#  109.04809    19.27074



# Test XeeX_serial -------------------------------------------------------------
# From ARE 212 Section 10: Non-Standard Standard Errors II
# and https://github.com/edrubin/ARE212

data(Wheat, package = "HistData")
wheatData <- Wheat %>% mutate(ones = 1) %>% na.omit() %>% as_tibble()

wheatModel <- lm(data = wheatData, Wheat ~ Wages)


# R Function
cannedNW <- sandwich::NeweyWest(wheatModel, lag = 3, prewhite = FALSE, sandwich = TRUE)
sqrt(diag(cannedNW))


# My Function
time = as.matrix(1:nrow(wheatData))
id = as.matrix(rep(1, nrow(wheatData)))
X = as.matrix(cbind(1,wheatData$Wages))
e = as.matrix(resid(wheatModel))

meat <- XeeX_serial(time = time, id = id, cutoff = 3, X = X, e = e, k=2)

sqrt(diag(solve(t(X) %*% X) %*% meat %*% solve(t(X) %*% X)))




# Test Conley-HAC SEs ----------------------------------------------------------
# https://darinchristensen.com/post/conley-correction/

df <- haven::read_dta("new_testspatial.dta") %>% drop_na(EmpClean00)

m <- fixest::feols(EmpClean00 ~ HDD + CDD | year + FIPS, data = df, demeaned = TRUE)

coords <- as.matrix(df[, c("latitude", "longitude")])
X <- m$X_demeaned
e <- m$residuals
time <- df$year
id <- df$FIPS
n <- length(e)
k <- ncol(X)
dist_cutoff <- 500
lag_cutoff <- 5


SE <- conley_ses(X=X, e=e, coords=coords, id=id, time=time, dist_cutoff=dist_cutoff, lag_cutoff=lag_cutoff)

sapply(SE, function(x) diag(sqrt(x))) %>% round(3)

## Should be this:
#     OLS     Spatial     Spatial_HAC
# HDD 0.650   0.886       0.721
# CDD 1.493   4.065       3.631



