library(tidyverse)

# setwd()
source("helper-conley.R")

# Test 
df <- haven::read_dta("new_testspatial.dta") %>% 
	filter(!is.na(EmpClean00))

m <- fixest::feols(EmpClean00 ~ HDD + CDD | year + FIPS, data = df, demeaned = TRUE)

M <- as.matrix(df[, c("latitude", "longitude")])
X <- m$X_demeaned
e <- m$residuals
time <- df$year
id <- df$FIPS
n <- length(e)
k <- ncol(M)


conley_ses(X, e, M, id, time, 100, 2, cores = 2)
