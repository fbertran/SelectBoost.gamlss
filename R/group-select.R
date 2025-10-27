
# Group Lasso (grpreg) selector for mu
select_grpreg_mu <- function(data, response, mu_scope_terms, penalty = c("grLasso","grMCP","grSCAD"), df_smooth = 6L) {
  penalty <- match.arg(penalty)
  if (!requireNamespace("grpreg", quietly = TRUE)) {
    stop("The 'grpreg' package is required for engine='grpreg'. Please install.packages('grpreg').")
  }
  gd <- .mu_group_design(data, mu_scope_terms, df_smooth = df_smooth)
  if (is.null(gd$X) || ncol(gd$X) == 0) return(character(0))
  X <- gd$X
  y <- data[[response]][gd$rows]
  keep <- which(!is.na(y) & is.finite(y))
  if (length(keep) == 0) return(character(0))
  if (length(keep) < length(y)) {
    X <- X[keep, , drop = FALSE]
    y <- y[keep]
  }
  if (nrow(X) == 0L) return(character(0))
  # Cross-validated group lasso
  cv <- grpreg::cv.grpreg(x = X, y = y, group = gd$groups, family = "gaussian",
                          penalty = penalty, seed = NULL)
  fit <- cv$fit
  beta <- stats::coef(fit, lambda = cv$lambda.min)[-1, , drop = FALSE]  # drop intercept
  # mark groups with any nonzero coefficient
  nonzero <- which(abs(beta[,1]) > 0)
  if (!length(nonzero)) return(character(0))
  grp_idx <- gd$groups[nonzero]
  sel_terms <- unique(gd$terms[grp_idx])
  sel_terms
}

# Sparse Group Lasso (SGL) selector for mu
select_sgl_mu <- function(data, response, mu_scope_terms, alpha = 0.95, df_smooth = 6L) {
  # alpha close to 1 emphasizes lasso within groups; alpha=0 is pure group lasso
  if (!requireNamespace("SGL", quietly = TRUE)) {
    stop("The 'SGL' package is required for engine='sgl'. Please install.packages('SGL').")
  }
  gd <- .mu_group_design(data, mu_scope_terms, df_smooth = df_smooth)
  if (is.null(gd$X) || ncol(gd$X) == 0) return(character(0))
  X <- gd$X
  y <- data[[response]][gd$rows]
  keep <- which(!is.na(y) & is.finite(y))
  if (length(keep) == 0) return(character(0))
  if (length(keep) < length(y)) {
    X <- X[keep, , drop = FALSE]
    y <- y[keep]
  }
  if (nrow(X) == 0L) return(character(0))
  # SGL needs a list with x and y; groups as integer vector starting at 1
  data_sgl <- list(x = X, y = as.numeric(y))
  idx <- gd$groups
  # Cross-validate
  cv <- SGL::cvSGL(data = data_sgl, index = idx, type = "linear", alpha = alpha, standardize = FALSE)
  fit <- cv$fit
  # pick lambda with min cv error
  lam <- cv$lam[which.min(cv$lldiff)]  # cv object stores cv curve; fallback to min index
  # coefficients at lam
  beta_mat <- SGL::SGL(data = data_sgl, index = idx, type = "linear", alpha = alpha, standardize = FALSE, lam = lam)$beta
  nonzero <- which(abs(beta_mat) > 0)
  if (!length(nonzero)) return(character(0))
  grp_idx <- idx[nonzero]
  sel_terms <- unique(gd$terms[grp_idx])
  sel_terms
}


# Heuristic link-like transform for positive / (0,1) / real responses
.linkish <- function(y) {
  if (all(is.finite(y)) && all(y > 0)) return(log(y))
  if (all(is.finite(y)) && all(y > 0) && all(y < 1)) return(qlogis(y))
  y
}


# Generic grpreg/SGL selectors for a parameter given a working response y_work
select_grpreg_param <- function(data, scope_terms, y_work, penalty = c("grLasso","grMCP","grSCAD"), df_smooth = 6L) {
  penalty <- match.arg(penalty)
  if (!requireNamespace("grpreg", quietly = TRUE)) stop("Package 'grpreg' is required for engine='grpreg'.")
  gd <- .param_group_design(data, scope_terms, df_smooth = df_smooth)
  if (is.null(gd$X) || ncol(gd$X) == 0) return(character(0))
  X <- gd$X
  yw <- as.numeric(y_work)[gd$rows]
  keep <- which(!is.na(yw) & is.finite(yw))
  if (length(keep) == 0) return(character(0))
  if (length(keep) < length(yw)) {
    X <- X[keep, , drop = FALSE]
    yw <- yw[keep]
  }
  if (nrow(X) == 0L) return(character(0))
  cv <- grpreg::cv.grpreg(x = X, y = yw, group = gd$groups, family = "gaussian",
                          penalty = penalty, seed = NULL)
  beta <- stats::coef(cv$fit, lambda = cv$lambda.min)[-1,,drop=FALSE]
  nz <- which(abs(beta[,1]) > 0)
  if (!length(nz)) return(character(0))
  sel <- unique(gd$terms[ gd$groups[nz] ])
  sel
}

select_sgl_param <- function(data, scope_terms, y_work, alpha = 0.95, df_smooth = 6L) {
  if (!requireNamespace("SGL", quietly = TRUE)) stop("Package 'SGL' is required for engine='sgl'.")
  gd <- .param_group_design(data, scope_terms, df_smooth = df_smooth)
  if (is.null(gd$X) || ncol(gd$X) == 0) return(character(0))
  X <- gd$X
  yw <- as.numeric(y_work)[gd$rows]
  keep <- which(!is.na(yw) & is.finite(yw))
  if (length(keep) == 0) return(character(0))
  if (length(keep) < length(yw)) {
    X <- X[keep, , drop = FALSE]
    yw <- yw[keep]
  }
  if (nrow(X) == 0L) return(character(0))
  dat <- list(x = X, y = as.numeric(yw))
  idx <- gd$groups
  cv <- SGL::cvSGL(data = dat, index = idx, type = "linear", alpha = alpha, standardize = FALSE)
  lam <- cv$lam[ which.min(cv$lldiff) ]
  fit <- SGL::SGL(data = dat, index = idx, type = "linear", alpha = alpha, standardize = FALSE, lam = lam)
  beta <- as.numeric(fit$beta)
  nz <- which(abs(beta) > 0)
  if (!length(nz)) return(character(0))
  sel <- unique(gd$terms[ idx[nz] ])
  sel
}

# Expand mu/sigma/nu/tau scope term labels into a model.matrix-ready RHS string.
# We replace common smoother constructors (pb(), cs()) with splines::bs() proxies for selection,
# so group lasso treats the whole basis as one group. Interactions with smoothers are not expanded.
.expand_terms_to_mm <- function(terms, df_smooth = 6L, converters = getOption("SelectBoost.gamlss.term_converters")) {
  if (length(terms) == 0L) return(character(0))
  if (!is.null(converters)) {
    if (!is.list(converters) || !all(vapply(converters, is.function, logical(1)))) {
      stop("`converters` must be a list of functions or NULL.")
    }
  }
  out <- character(length(terms))
  for (i in seq_along(terms)) {
    t <- terms[i]
    # pb(x) -> splines::bs(x, df = df_smooth)
    t <- gsub("(\\bpb\\(\\s*([[:alnum:]_\\.]+)\\s*\\))",
              sprintf("splines::bs(\\2, df=%d)", df_smooth), t, perl = TRUE)
    # cs(x) -> splines::bs(x, df = df_smooth)
    t <- gsub("(\\bcs\\(\\s*([[:alnum:]_\\.]+)\\s*\\))",
              sprintf("splines::bs(\\2, df=%d)", df_smooth), t, perl = TRUE)
    # pbm(x) -> splines::bs(x, df = df_smooth)
    t <- gsub("(\\bpbm\\(\\s*([[:alnum:]_\\.]+)\\s*\\))",
              sprintf("splines::bs(\\2, df=%d)", df_smooth), t, perl = TRUE)
    # lo(x, ...) -> splines::ns(x, df = df_smooth) using first argument
    t <- gsub("(\\blo\\(\\s*([[:alnum:]_\\.]+)[^)]*\\))",
              sprintf("splines::ns(\\2, df=%d)", df_smooth), t, perl = TRUE)
    if (!is.null(converters) && length(converters)) {
      for (conv in converters) {
        res <- conv(t, df_smooth = df_smooth)
        if (!is.null(res)) {
          if (!is.character(res) || length(res) != 1L) {
            stop("Custom term converters must return a single character string or NULL.")
          }
          t <- res
        }
      }
    }
    out[i] <- t
  }
  out
}

.mu_group_design <- function(data, mu_scope_terms, df_smooth = 6L) {
  mm_terms <- .expand_terms_to_mm(mu_scope_terms, df_smooth = df_smooth)
  f <- stats::as.formula(paste("~ 0 +", paste(mm_terms, collapse = " + ")))
  X <- stats::model.matrix(f, data)
  
  tl <- mu_scope_terms
  
  rows <- rownames(X)
  if (is.null(rows)) {
    row_idx <- seq_len(nrow(X))
  } else if (!is.null(rownames(data))) {
    row_idx <- match(rows, rownames(data))
  } else {
    row_idx <- suppressWarnings(as.integer(rows))
  }
  if (anyNA(row_idx) || length(row_idx) != nrow(X)) {
    row_idx <- seq_len(nrow(X))
  }
  row_idx <- as.integer(row_idx)
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  rownames(X) <- NULL
  # groups: one id per original term, repeated for its columns
  # build using attr(terms, "term.labels") mapping

  assign <- attr(X, "assign")
  if (is.null(assign)) {
    assign <- integer(ncol(X))
  }
  groups <- assign
  
  list(X = X, groups = groups, terms = tl, rows = row_idx)
}

.param_group_design <- function(data, scope_terms, df_smooth = 6L) {
  mm_terms <- .expand_terms_to_mm(scope_terms, df_smooth = df_smooth)
  f <- stats::as.formula(paste("~ 0 +", paste(mm_terms, collapse = " + ")))
  X <- stats::model.matrix(f, data)
  
  tl <- scope_terms
  
  rows <- rownames(X)
  if (is.null(rows)) {
    row_idx <- seq_len(nrow(X))
  } else if (!is.null(rownames(data))) {
    row_idx <- match(rows, rownames(data))
  } else {
    row_idx <- suppressWarnings(as.integer(rows))
  }
  if (anyNA(row_idx) || length(row_idx) != nrow(X)) {
    row_idx <- seq_len(nrow(X))
  }
  row_idx <- as.integer(row_idx)

  X <- as.matrix(X)
  storage.mode(X) <- "double"
  rownames(X) <- NULL
  
  assign <- attr(X, "assign")
  if (is.null(assign)) {
    assign <- integer(ncol(X))
  }
  groups <- assign
  
  list(X = X, groups = groups, terms = tl, rows = row_idx)
}
