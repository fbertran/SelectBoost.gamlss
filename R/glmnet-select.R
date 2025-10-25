
# Experimental glmnet-based selection for the mu parameter
# Requires: glmnet
# Limitations: currently supports predictors that can be expanded via model.matrix.

.glmnet_prepare_response <- function(y, family) {
  if (is.null(y)) {
    stop("glmnet selectors require a non-null response vector.")
  }
  if (identical(family, "binomial")) {
    if (is.factor(y)) {
      lev <- levels(y)
      if (length(lev) != 2L) {
        stop("binomial glmnet selector requires a two-level factor response.")
      }
      return(as.numeric(y == lev[length(lev)]))
    }
    if (is.logical(y)) {
      return(as.integer(y))
    }
    y_num <- as.numeric(y)
    uniq <- sort(unique(y_num))
    uniq <- uniq[is.finite(uniq)]
    if (!length(uniq)) {
      stop("binomial glmnet selector received a response with no finite values.")
    }
    if (length(uniq) == 1L) {
      if (!uniq %in% c(0, 1)) {
        stop("binomial glmnet selector requires values coded as 0/1.")
      }
      return(y_num)
    }
    if (length(uniq) == 2L) {
      if (all(uniq %in% c(0, 1))) {
        return(y_num)
      }
      # map the smaller value to 0 and the larger to 1
      return(as.numeric(y_num == max(uniq)))
    }
    stop("binomial glmnet selector requires a binary response (0/1 or two-level factor).")
  }
  as.numeric(y)
}

.glmnet_select_terms <- function(data, response, mu_scope_terms, alpha = 1, family = c("gaussian","binomial","poisson")) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("glmnet is required for engine='glmnet'. Please install.packages('glmnet').")
  }
  family <- match.arg(family)
  # Keep numeric columns that appear as simple terms
  vars <- intersect(mu_scope_terms, names(data))
  if (!length(vars)) return(character(0))

  term_labels <- sprintf("`%s`", vars)
  mm_formula <- stats::as.formula(paste("~ 0 +", paste(term_labels, collapse = " + ")))
  X <- stats::model.matrix(mm_formula, data = data)
  y <- .glmnet_prepare_response(data[[response]], family)
  
  cv <- glmnet::cv.glmnet(x = X, y = y, alpha = alpha, family = family, standardize = FALSE)
  fit <- glmnet::glmnet(x = X, y = y, alpha = alpha, family = family, lambda = cv$lambda.min, standardize = FALSE)
  beta <- as.matrix(coef(fit))[-1, , drop = FALSE]  # drop intercept
  sel <- rownames(beta)[abs(beta[,1]) > 0]
  if (!length(sel)) return(character(0))
  assign_vec <- attr(X, "assign")
  if (is.null(assign_vec)) {
    base_vars <- vars
  } else {
    nz <- match(sel, colnames(X), nomatch = 0L)
    nz <- nz[nz > 0L]
    if (!length(nz)) return(character(0))
    term_idx <- assign_vec[nz]
    term_idx <- term_idx[term_idx >= 0 & term_idx <= length(vars)]
    base_vars <- unique(vars[term_idx])
  }
  intersect(base_vars, vars)
}

# convenience wrappers
select_lasso_mu  <- function(data, response, mu_scope_terms, family = c("gaussian","binomial","poisson"))
  .glmnet_select_terms(data, response, mu_scope_terms, alpha = 1, family = family)
select_ridge_mu  <- function(data, response, mu_scope_terms, family = c("gaussian","binomial","poisson"))
  .glmnet_select_terms(data, response, mu_scope_terms, alpha = 0, family = family)
select_glmnet_mu <- function(data, response, mu_scope_terms, alpha = 0.5, family = c("gaussian","binomial","poisson"))
  .glmnet_select_terms(data, response, mu_scope_terms, alpha = alpha, family = family)
