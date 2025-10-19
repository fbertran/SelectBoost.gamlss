
# Experimental glmnet-based selection for the mu parameter
# Requires: glmnet
# Limitations: currently supports numeric predictors in mu_scope, Gaussian family approximation.

.glmnet_select_terms <- function(data, response, mu_scope_terms, alpha = 1) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("glmnet is required for engine='glmnet'. Please install.packages('glmnet').")
  }
  # Keep numeric columns that appear as simple terms
  vars <- intersect(mu_scope_terms, names(data))
  vars <- vars[vapply(data[vars], is.numeric, logical(1))]
  if (!length(vars)) return(character(0))

  X <- stats::model.matrix(stats::as.formula(paste("~ 0 +", paste(vars, collapse = "+"))), data = data)
  y <- data[[response]]
  # gaussian path
  cv <- glmnet::cv.glmnet(x = X, y = y, alpha = alpha, family = "gaussian", standardize = FALSE)
  fit <- glmnet::glmnet(x = X, y = y, alpha = alpha, family = "gaussian", lambda = cv$lambda.min, standardize = FALSE)
  beta <- as.matrix(coef(fit))[-1, , drop = FALSE]  # drop intercept
  sel <- rownames(beta)[abs(beta[,1]) > 0]
  # map back to base variable names (drop factor levels if any)
  base_vars <- unique(sub("(.*?)([:].*)?$", "\\1", sel))
  intersect(base_vars, vars)
}

# convenience wrappers
select_lasso_mu  <- function(data, response, mu_scope_terms) .glmnet_select_terms(data, response, mu_scope_terms, alpha = 1)
select_ridge_mu  <- function(data, response, mu_scope_terms) .glmnet_select_terms(data, response, mu_scope_terms, alpha = 0)
select_glmnet_mu <- function(data, response, mu_scope_terms, alpha = 0.5) .glmnet_select_terms(data, response, mu_scope_terms, alpha = alpha)
