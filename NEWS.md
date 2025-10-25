# SelectBoost.gamlss 0.2.0

### Highlights
- **Grouped selection for all parameters (μ, σ, ν, τ)** with `engine = "grpreg"` (group lasso/MCP/SCAD) and `engine = "sgl"` (sparse group lasso), including **factors, splines (`pb()`/`cs()`), and interactions** treated as *single groups*.
- **Per-parameter engines**: `engine`, `engine_sigma`, `engine_nu`, `engine_tau` can be mixed (stepwise / glmnet / grpreg / sgl).
- **Glmnet support** (lasso/ridge/elastic-net) extended beyond μ via working-response proxies for σ/ν/τ.
- **Glmnet selectors** now accept `glmnet_family` (gaussian/binomial/poisson) and handle factor predictors via model-matrix expansion.
- **Tuning framework**: `tune_sb_gamlss()` with **stability** or **deviance** metrics (K-fold), **progress bars**, and a **complexity penalty**.
- **Fast deviance paths** for common families (auto-used in deviance CV): `NO`, `PO`, `LOGNO`, `GA`, `IG`, `LO`, `LOGITNO`, `GEOM`, `BE`, `NBI`, `NBII`, `BI`, and native shortcuts via `gamlss.dist` for many others (e.g., `LOGLOG`, `DEL`, `ZAGA`, `ZIP/ZIP2`, `ZAIG`, `ZALG`, `ZIBI/ZIBB`, `PARETO`, `SEP1/SEP2`, `ZIPF/ZIPFmu`, `BCT`, `BCPE`, `SICHEL`, `GLG`, `BETA4`, `RS`, `WEI`, `GIG`), with graceful fallbacks.
- **Fast deviance** dynamically calls `gamlss.dist::d<family>()` when available, broadening zero-inflated/hurdle coverage without manual whitelists.
- **Group knockoffs** (approximate) for FDR-style control: `knockoff_filter_mu()`, `knockoff_filter_param()`.
- **Compiled speedups** (Rcpp/RcppArmadillo) for scaling/cor; **parallel bootstraps** via `future.apply`.

### New user-facing functions
- `tune_sb_gamlss()`, `knockoff_filter_mu()`, `knockoff_filter_param()`
- `fast_vs_generic_ll()`, `check_fast_vs_generic()`
- `effect_plot()` (quick partial effect visualizer for the final selected model)

### New arguments in `sb_gamlss()`
- `engine_sigma`, `engine_nu`, `engine_tau` — choose engines per-parameter
- `grpreg_penalty` (`grLasso`/`grMCP`/`grSCAD`), `sgl_alpha`
- `df_smooth` — basis size for grouped-smoother proxies
- `progress` — progress bar for sequential bootstraps
- (still) `glmnet_alpha` — 0=ridge, 1=lasso, (0,1)=EN
- `glmnet_family` — choose gaussian/binomial/poisson for glmnet selectors

### Documentation & vignettes
- **Real Data Examples** (including **growth/BCT** on `gamlss.data::boys`)
- **Advanced Real Data Examples** (ZIP/ZINB on `bioChemists`, ZAGA on `airquality::Ozone`, longitudinal growth on `nlme::Orthodont` with random intercepts)
- **Benchmarks** (engine timings) and **Fast deviance microbenchmarks**
- **Fast deviance equality** (accuracy checks) + **wide-family sweep** with per-family tolerances & skip reasons
- **pkgdown** site scaffolding + GitHub Actions workflow
- README quick start now covers factor effects, BCT four-parameter example, and deviance-based tuning metrics.

### Testing & quality
- Unit tests for fast vs generic deviance (accuracy and presence of native densities)
- Opt-in **long tests** via `options(SelectBoost.gamlss.run_long_tests=TRUE)` or `RUN_LONG_TESTS=true`

### Notes
- Some grouped/knockoff features are optional and require packages in **Suggests** (`grpreg`, `SGL`, `knockoff`, `glmnet`, etc.).
- Smooths are proxied with `splines::bs(df = df_smooth)` *for selection only*; the final `gamlss` fit remains as specified.

# SelectBoost.gamlss 0.1.0

- First draft: bootstrap stability-selection over GAMLSS parameters (mu/sigma/nu/tau).
- Optional pre-standardization of numeric predictors (stored for prediction).
- AICc helper.
- Plotting + prediction helpers.
