# Internal helpers bridging gamlss design matrices with SelectBoost machinery

# Build a model matrix for a RHS-only formula, sanitize column names, and
# return mappings back to original term labels.
.sb_scope_matrix <- function(data, scope_formula, prefix) {
  if (is.null(scope_formula)) {
    return(list(
      matrix = NULL,
      original_colnames = character(0),
      sanitized_colnames = character(0),
      term_map = character(0),
      term_labels = character(0),
      assign = integer(0)
    ))
  }
  if (!inherits(scope_formula, "formula")) {
    stop("scope_formula must be a formula")
  }
  # Always drop the intercept; SelectBoost works on centred columns
  scope_no_intercept <- try(stats::update(scope_formula, ~ . - 1), silent = TRUE)
  if (inherits(scope_no_intercept, "try-error")) scope_no_intercept <- scope_formula
  mm <- stats::model.matrix(scope_no_intercept, data)
  if (!is.matrix(mm) || !ncol(mm)) {
    return(list(
      matrix = NULL,
      original_colnames = character(0),
      sanitized_colnames = character(0),
      term_map = character(0),
      term_labels = attr(stats::terms(scope_no_intercept), "term.labels"),
      assign = integer(0)
    ))
  }
  original_colnames <- colnames(mm)
  assign <- attr(mm, "assign")
  term_labels <- attr(stats::terms(scope_no_intercept), "term.labels")
  # Sanitise column names, keeping a deterministic prefix for debugging
  base_names <- make.names(original_colnames, unique = FALSE)
  sanitized <- make.unique(paste0(prefix, "__", base_names), sep = "__")
  colnames(mm) <- sanitized
  term_map <- if (length(assign)) {
    tl <- if (length(term_labels)) term_labels else character(0)
    out <- rep("", length(assign))
    nz <- assign > 0 & assign <= length(tl)
    out[nz] <- tl[assign[nz]]
    out
  } else character(0)
  if (length(term_map)) names(term_map) <- sanitized
  list(
    matrix = mm,
    original_colnames = original_colnames,
    sanitized_colnames = sanitized,
    term_map = term_map,
    term_labels = term_labels,
    assign = assign
  )
}

# Prepare SelectBoost correlated-resampling artefacts for a scope model matrix.
.sb_prepare_selectboost <- function(
    data,
    scope_formula,
    prefix,
    B,
    corr_func = "cor",
    group_fun = SelectBoost::group_func_2,
    corr_threshold = 0.5,
    use_groups = TRUE
) {
  scope_mm <- .sb_scope_matrix(data, scope_formula, prefix)
  X <- scope_mm$matrix
  if (is.null(X) || !ncol(X)) {
    empty <- matrix(0, nrow = NROW(data), ncol = 0L)
    colnames(empty) <- character(0)
    attr(empty, "nosimul") <- TRUE
    attr(empty, "colstosimul") <- integer(0)
    attr(empty, "nsimul") <- B
    return(list(
      matrix = empty,
      simulations = empty,
      colstosimul = integer(0),
      groups = NULL,
      corr = NULL,
      corr_sign = NULL,
      adjust = NULL,
      scope = scope_mm
    ))
  }
  Xnorm <- SelectBoost::boost.normalize(X)
  colnames(Xnorm) <- scope_mm$sanitized_colnames
  corr <- SelectBoost::boost.compcorrs(Xnorm, corrfunc = corr_func)
  corr_sign <- SelectBoost::boost.correlation_sign(corr)
  groups <- if (isTRUE(use_groups)) {
    SelectBoost::boost.findgroups(corr, group = group_fun, corr = corr_threshold)
  } else {
    g_groups <- lapply(seq_len(ncol(Xnorm)), function(j) j)
    attr(g_groups, "type") <- "singleton"
    attr(g_groups, "length.groups") <- rep(1L, length(g_groups))
    list(groups = g_groups)
  }
  adjust <- SelectBoost::boost.adjust(Xnorm, groups$groups, corr_sign)
  sims <- SelectBoost::boost.random(Xnorm, adjust$Xpass, adjust$vmf.params, B = B)
  if (is.null(attr(sims, "colstosimul"))) {
    attr(sims, "colstosimul") <- integer(0)
  }
  list(
    matrix = Xnorm,
    simulations = sims,
    colstosimul = attr(sims, "colstosimul"),
    groups = groups,
    corr = corr,
    corr_sign = corr_sign,
    adjust = adjust,
    scope = scope_mm
  )
}
