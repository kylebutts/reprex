## -----------------------------------------------------------------------------
## oaxaca-blinder-estimator.R
## Kyle Butts, CU Boulder Economics 
## 
## This creates R functions for an Oaxaca-Blinder Estimator following:
## Patrick Kline - Oaxaca-Blinder as a Reweighting Estimator (2011)
## Patrick Kline - A note on variance estimation for the Oaxaca estimator of average treatment effects (2014)
## -----------------------------------------------------------------------------


## Kline (2011) ----------------------------------------------------------------

#' @description This function returns the Oaxaca-Blinder Estimand for the Average Treatment Effect on the Treated, following Patrick Kline - Oaxaca-Blinder as a Reweighting Estimator (2011)
#'
#' @param data dataframe to estimate with
#' @param formula Formula for linear regression specification. Either formula or string object. Do not include the treatment variable in this
#' @param treat string for treatment variable, must be 1 = treat, 0 = control
#' @param cluster optional - string for cluster variable. Variable "indicating which observations belong to the same cluster". Passed along to `clubSandwich::vcovCR()`
oaxaca_estimate <- function(data, formula, treat, cluster = NA) {
	
	## Regression using Control units ------------------------------------------
	
	# Prepare formula
	if(is.character(formula)) formula <- as.formula(formula)
	Y_var <- all.vars(formula[[2]])
	X_vars <- all.vars(formula)[-1]
	
	# Subset by control and treated
	data_control <- data[data[[treat]] == 0, ]
	data_treat <- data[data[[treat]] == 1, ]
	control_reg <- lm(formula, data = data_control) 
	
	
	## Mean of TE --------------------------------------------------------------
	
	# Counterfactual Mean
	pred <- predict(control_reg, data)
	
	# Y_i(1) - \hat{Y}_i(0)
	diff <- data[[Y_var]] - pred
	
	# E[Y_i(1) - \hat{Y}_i(0) | D_i = 1]
	te <- mean(diff[data[[treat]] == 1])
	
	
	## Variance of TE ----------------------------------------------------------
	
	# Stata Robust Var-Cov matrix
	if(is.na(cluster)) {
		V0 <- sandwich::vcovHC(control_reg, type = "HC1")
	} else {
		V0 <- clubSandwich::vcovCR(control_reg, cluster = data_control[[cluster]], type = "CR1S")
	}
	
	# Var(te) comes from out of sample prediction variance formula
	X <- model.matrix(formula, data = data)
	D <- data[[treat]]
	V1 <- var(data_treat[[Y_var]]) / sum(D)
	Vdiff <- V1 + (t(D) %*% X %*% V0 %*% t(X) %*% D)/(sum(D)^2)
	
	se <- as.vector(sqrt(Vdiff))
	
	return(list("te" = te, "se" = se))	
}




## Kline (2014) ----------------------------------------------------------------

#' @description This function returns the Oaxaca-Blinder Estimand for the Average Treatment Effect on the Treated, following Patrick Kline - A note on variance estimation for the Oaxaca estimator of average treatment effects (2014). This is a more robust version of the estimator. 
#'
#' @param data dataframe to estimate with
#' @param formula Formula for linear regression specification. Either formula or string object. Do not include the treatment variable in this
#' @param treat string for treatment variable, must be 1 = treat, 0 = control
#' @param cluster optional - string for cluster variable. Variable "indicating which observations belong to the same cluster". Passed along to `clubSandwich::vcovCR()`
oaxaca_estimate_robust <- function(data, formula, treat, cluster = NA) {
	
	## Regression using Control units ------------------------------------------
	
	# Prepare formula
	if(is.character(formula)) formula <- as.formula(formula)
	Y_var <- all.vars(formula[[2]])
	X_vars <- all.vars(formula)[-1]
	
	D <- data[[treat]]
	X <- model.matrix(formula, data = data)
	
	
	
	## Step 1. (Kline 2014) ----------------------------------------------------
	# Subset by control and treated
	data_control <- data[data[[treat]] == 0, ]
	data_treat <- data[D == 1, ]
	control_reg <- lm(formula, data = data_control) 
	
	
	## Step 2. (Kline 2014) ----------------------------------------------------
	pred <- predict(control_reg, data)
	Y_star <- D * (data[[Y_var]] - pred) + (1-D) * data[[Y_var]]
	
	## Step 3. (Kline 2014) ----------------------------------------------------
	data_update <- as.data.frame(cbind(D, (1-D) * X))
	data_update[[Y_var]] <- Y_star
	
	update_formula <- update(formula, paste("~ . +", "D - 1"))
	
	reg_step_3 <- lm(update_formula, data = data_update)
	
	
	## Step 4. (Kline 2014) ----------------------------------------------------
	te <- coef(reg_step_3)[["D"]]
	
	# Stata Robust Var-Cov matrix
	if(is.na(cluster)) {
		V <- sandwich::vcovHC(reg_step_3, type = "HC1")
	} else {
		V <- clubSandwich::vcovCR(
				reg_step_3, 
				cluster = data[[cluster]], 
				type = "CR1S"
			)
	}
	
	V_theta_hat <- V["D","D"]
	V_beta <- V[setdiff(rownames(V),"D"),setdiff(rownames(V),"D")]
	V_beta_theta <- V["D",setdiff(rownames(V),"D")]
	
	
	## Step 5. (Kline 2014) ----------------------------------------------------
	
	mu_X_treated <- colMeans(data[D == 1, X_vars])
	
	V_theta <- V_theta_hat + 
		t(mu_X_treated) %*% V_beta %*% mu_X_treated -
		2 * t(mu_X_treated) %*% V_beta_theta
	
	se <- as.vector(sqrt(V_theta))
	
	
	## Return estimate and se --------------------------------------------------
	return(list("te" = te, "se" = se))	
}

