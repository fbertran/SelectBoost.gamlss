
# Knockoff-based (approximate) group FDR control for mu/parameters
# This uses variable-level knockoffs and aggregates to groups by max |W|.

#' Knockoff filter for mu (approximate group control)
#' @param data data.frame
#' @param response response variable name
#' @param mu_scope RHS-only term labels
#' @param fdr target FDR level
#' @param df_smooth df for smoother proxies (splines::bs)
#' @return character vector of selected term names
#' @export
knockoff_filter_mu <- function(data, response, mu_scope, fdr = 0.1, df_smooth = 6L) {
  if (!requireNamespace("knockoff", quietly = TRUE)) stop("Package 'knockoff' required for knockoff filtering.")
  gd <- .mu_group_design(data, attr(stats::terms(mu_scope), "term.labels"), df_smooth = df_smooth)
  if (is.null(gd$X) || ncol(gd$X) == 0) return(character(0))
  y <- data[[response]]
  K <- knockoff::create.fixed(gd$X)
  W <- knockoff::stat.glmnet_coefdiff(gd$X, K, y)
  t <- knockoff::knockoff.threshold(W, fdr = fdr, offset = 1)
    # group-level statistics: W_g = max_j |W_j| within group, signed by argmax
  grp_ids <- sort(unique(gd$groups))
  Wg <- sapply(grp_ids, function(g) {
    idx <- which(gd$groups == g)
    jj <- idx[which.max(abs(W[idx]))]
    sign(W[jj]) * max(abs(W[idx]))
  })
  names(Wg) <- grp_ids
  tg <- knockoff::knockoff.threshold(Wg, fdr = fdr, offset = 1)
  keep <- which(Wg >= tg)
  if (!length(keep)) return(character(0))
  sel_groups <- grp_ids[keep]
  unique(gd$terms[match(sel_groups, gd$groups)])
}

#' Knockoff filter for sigma/nu/tau (approximate group control)
#' @param data data.frame
#' @param scope RHS-only term labels
#' @param y_work working response (numeric)
#' @param fdr target FDR level
#' @param df_smooth df for smoother proxies
#' @export
knockoff_filter_param <- function(data, scope, y_work, fdr = 0.1, df_smooth = 6L) {
  if (!requireNamespace("knockoff", quietly = TRUE)) stop("Package 'knockoff' required for knockoff filtering.")
  gd <- .param_group_design(data, attr(stats::terms(scope), "term.labels"), df_smooth = df_smooth)
  if (is.null(gd$X) || ncol(gd$X) == 0) return(character(0))
  K <- knockoff::create.fixed(gd$X)
  W <- knockoff::stat.glmnet_coefdiff(gd$X, K, y_work)
  t <- knockoff::knockoff.threshold(W, fdr = fdr, offset = 1)
    # group-level statistics: W_g = max_j |W_j| within group, signed by argmax
  grp_ids <- sort(unique(gd$groups))
  Wg <- sapply(grp_ids, function(g) {
    idx <- which(gd$groups == g)
    jj <- idx[which.max(abs(W[idx]))]
    sign(W[jj]) * max(abs(W[idx]))
  })
  names(Wg) <- grp_ids
  tg <- knockoff::knockoff.threshold(Wg, fdr = fdr, offset = 1)
  keep <- which(Wg >= tg)
  if (!length(keep)) return(character(0))
  sel_groups <- grp_ids[keep]
  unique(gd$terms[match(sel_groups, gd$groups)])
}
