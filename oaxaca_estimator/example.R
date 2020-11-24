## Example from Kline (2011) using Lalonde Data --------------------------------

# setwd()

# Load functions
source("helper-oaxaca.R")

# Load Data
data <- haven::read_dta("cps3re74.dta")

formula <- "re78 ~ 1 + age + age2 + ed + black + hisp + married + nodeg + re75 + re74"

# Oaxaca-Blindor Estimator
oaxaca_estimate(data, formula, "treat")
oaxaca_estimate_robust(data, formula, "treat")

# Cluster by education
oaxaca_estimate(data, formula, "treat", "ed")
oaxaca_estimate_robust(data, formula, "treat", "ed")

