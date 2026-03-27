# Model Descriptions

Narrative descriptions of each model used in:
> Gabauer, Gupta, Karmakar & Nielsen (2022). *Stock Market Bubbles and the Forecastability of Gold Returns (and Volatility)*. University of Pretoria Working Paper 2022-28.

For the mathematical equations see [`model_spec.md`](model_spec.md).

---

## 1. LPPLS — Log-Periodic Power Law Singularity Model

**Role:** The fundamental building block. Fitted to a window of the log price-dividend ratio of a single stock market to determine whether a bubble pattern is present and, if so, when it is expected to terminate.

**What it does:** The LPPLS model captures the signature of a speculative bubble: superexponential price growth (or decay for negative bubbles) decorated with log-periodic oscillations that increase in frequency as the system approaches a critical point. When fitted to a time window of a stock index's log price-dividend ratio, the model returns seven parameters. The most diagnostic is `B`: a negative `B` means the series accelerates upward toward a crash (positive/crash bubble), while a positive `B` means it decelerates downward toward a rally (negative/recovery bubble). The critical time `tc` marks the expected crash or rally date.

**Why it is used instead of alternatives:** Standard bubble detection tests (e.g., PSY, SADF) detect departure from a fundamental value but cannot characterise the *timing* or *multi-scale structure* of a bubble. The LPPLS framework, grounded in the Johansen-Ledoit-Sornette (JLS) model of rational expectation bubbles with crash hazard rates, produces both a bubble indicator and a termination forecast. The Filimonov-Sornette (2013) reparameterisation used here eliminates the phase parameter `φ`, separating the problem into an analytic 4×4 linear solve and a lower-dimensional nonlinear search, which is numerically more stable and faster than fitting the original 7-parameter form directly.

**Input:** A time window `[t1, t2]` of the log price-dividend ratio for one country (daily observations, N points).

**Output:** Seven parameter estimates `{tc, m, ω, A, B, C1, C2}` and a fit quality score. The fit is either accepted or rejected by the filter conditions (see model_spec.md §2.3).

---

## 2. Multi-Scale LPPLS Confidence Indicator (LPPLS-CI)

**Role:** Aggregates hundreds of individual LPPLS fits into a single number per time point per scale, quantifying how strongly a bubble pattern is present in the data right now. This is the output that feeds into the forecasting model as predictors.

**What it does:** Rather than fitting LPPLS to a single fixed window, the CI sweeps over many overlapping windows of varying length ending at the current date `t2`. For each window `[t1, t2]` the LPPLS model is calibrated and the fit is "qualified" only if it passes all five filter conditions (bounds on `m`, `ω`, `tc`, oscillation count `O`, and damping ratio `D`). The CI is the fraction of qualified fits. This is done separately for three window-length ranges (short, medium, long-term), yielding three scales. Each scale is also split into positive and negative bubble indicators. The result is 6 CI values per country per week.

**Why three scales:** Different classes of market participants react to information at different horizons — traders to short-term signals, investors to long-term ones (Heterogeneous Market Hypothesis, Müller et al. 1997). Short-term CI values tend to precede medium-term ones, which precede long-term ones, as a bubble matures across scales. Using all three scales independently preserves this lead-lag information.

**Why separate positive and negative:** Crashes and recoveries carry different information for gold as a safe haven. Positive bubble indicators (crash precursors) tend to be stronger predictors of gold returns than negative ones, because investors flee to gold primarily during equity market turmoil.

**Input:** Full history of the log price-dividend ratio for one country up to the current date.

**Output:** 6 weekly time series per country (positive/negative × short/medium/long). Across 7 G7 countries this gives 42 predictor series; across all 12 countries (G7 + BRICS) it gives 72.

---

## 3. AR(1) Benchmark — Model 1

**Role:** The null model against which all forecasting gains are measured. Represents the "no bubble information" baseline.

**What it does:** A simple first-order autoregression fits the current gold log-return as a linear function of last week's return plus an intercept. It has only two parameters (`a0`, `a1`) and captures any short-term momentum or mean-reversion in gold returns.

**Why it is the benchmark:** It is the simplest possible time-series model that respects the autocorrelation structure of returns. Using a more complex benchmark (e.g., AR(p), ARMA) would make it harder to attribute any forecasting gain to the LPPLS-CI indicators. Results for higher-order AR models are noted in the paper to be qualitatively similar. The Clark–West test is specifically designed for nested model comparison, making the AR(1) the natural benchmark.

**Input:** Weekly gold log-return series.

**Output:** One-step-ahead (and h-step block-averaged) forecasts of weekly gold log-returns.

---

## 4. LASSO + Block-Mean Forecasting — Model 2

**Role:** The main forecasting model. Incorporates the LPPLS-CI predictors to produce point forecasts of gold returns (and volatility) that beat the AR(1) benchmark at short-to-medium horizons.

**What it does:** An AR(1) model is augmented with up to 72 LPPLS-CI predictors. Because the number of predictors exceeds what ordinary least squares can handle reliably, LASSO (Least Absolute Shrinkage and Selection Operator) is used. LASSO imposes an L1 penalty on all coefficients (including `a0` and `a1`), automatically shrinking irrelevant predictors to zero and selecting the informative ones. The regularisation strength `λ` is chosen by cross-validation (`cv.glmnet`, `lambda.min`).

**Why LASSO and not OLS or ridge:** Gold returns have heavy upper tails (Jarque-Bera statistic of 4878 in the paper), which means parametric forecasting methods fail. LASSO's hard zeroing of predictors performs well in high-dimensional, heavy-tailed settings. Ridge regression would retain all predictors; Lasso's sparsity is more interpretable and empirically superior here according to Karmakar et al. (2021a).

**Block-mean correction:** Standard LASSO residuals are serially correlated. A quantile-based block-mean correction (Karmakar et al. 2021a, extending Zhou et al. 2010) averages the residuals over blocks of length `h` and adds this mean back to the LASSO fitted values. This produces point forecasts of the `h`-week-ahead mean return with guaranteed accuracy properties under the law of large numbers, even for heavy-tailed data, as long as `h` is moderately large.

**Evaluation — Clark–West (CW) test:** Because Model 2 nests Model 1, standard forecast comparison tests (Diebold–Mariano) are invalid. The Clark–West test adjusts for the finite-sample bias introduced by the larger model and produces valid p-values for the one-sided test that Model 2 beats Model 1.

**Evaluation — Diebold–Mariano (DM) test:** Used only to compare the positive-only vs. negative-only versions of Model 2 against each other (non-nested comparison).

**Input:** Weekly gold log-returns + LPPLS-CI predictor matrix (42 or 72 columns).

**Output:** Rolling pseudo-out-of-sample point forecasts at horizons `h ∈ {1, 2, 4, 8, 12}` weeks, for training lengths `T ∈ {100, 250, 500, 1000}`.

---

## 5. TVPGARCH — Time-Varying Parameter GARCH(1,1)

**Role:** Produces the volatility series that serves as the dependent variable when forecasting gold return *volatility* (rather than returns). It is not itself a forecasting model — it is a pre-processing step that extracts a time-varying volatility estimate from the gold return series.

**What it does:** A standard GARCH(1,1) model assumes constant parameters across the full sample. For a nearly 50-year weekly series of gold returns, this is implausible — the relationship between past squared returns and today's variance clearly shifts over time (e.g., the 1980 gold price spike, the 2008 GFC). The TVPGARCH model allows the three GARCH parameters (intercept `α0`, ARCH coefficient `α1`, GARCH coefficient `β1`) to be smooth functions of rescaled time `t/n` rather than fixed constants. These functions are estimated by local kernel-weighted maximum likelihood, where observations closer in time to the target date receive higher weight (Epanechnikov kernel, bandwidth `b_n = 0.25`). The resulting `σ²_t` series is then used as the target variable in place of raw gold returns when fitting Models 1 and 2 for volatility forecasting.

**Why time-varying parameters:** The gold volatility series shown in Figure 1 of the paper exhibits clear non-stationarity — it spikes to ~160 in the early 1980s and stays near zero for long stretches. A constant-parameter GARCH would be severely misspecified. The TVPGARCH of Karmakar and Roy (2021) has been shown to outperform constant-parameter GARCH in this setting.

**Why Epanechnikov kernel:** The Epanechnikov kernel is the optimal kernel in the mean-squared-error sense for local polynomial smoothing. It assigns a parabolic weight profile to nearby observations and exact zero weight beyond the bandwidth, which is computationally efficient.

**Input:** Weekly gold log-return series.

**Output:** Weekly conditional variance series `σ²_t` for 1973–2020.

---

## 6. TVP-VAR / Total Connectedness Index (TCI)

**Role:** An optional alternative aggregation of the 6×7 = 42 (or 6×12 = 72) individual LPPLS-CI series into 6 compact connectedness measures — one per bubble indicator type (pos/neg × short/medium/long). The TCI series are used in place of (or alongside) the raw CI indicators as predictors in Model 2, to test whether network-level bubble interconnectedness carries forecasting power beyond individual country indicators.

**What it does:** A time-varying parameter vector autoregression (TVP-VAR) is fitted to the 7 (or 12) LPPLS-CI values of one bubble type at one scale at each point in time. The TVP-VAR has time-varying coefficient matrices and covariance matrices, estimated via the Kalman filter with random-walk state evolution. The model is then transformed to its VMA representation, from which generalised forecast error variance decompositions (GFEVD) are computed. The Total Connectedness Index (TCI) summarises how much of the forecast error variance of each country's bubble indicator is explained by shocks to other countries' indicators, averaged over all countries and normalised to [0,1]. A high TCI means bubble signals are highly synchronised across countries (systemic event); a low TCI means idiosyncratic national signals.

**Why connectedness matters:** The paper shows that LPPLS-CI values across G7 and BRICS countries are highly synchronised around major global events (Black Monday, Dot-com, GFC, 2015-16 sell-off). If this synchronisation itself carries information about gold's safe-haven demand, the TCI — which captures the degree of cross-country bubble spillover — may be an informative predictor even when individual country indicators are uninformative. Empirically (Table 2), TCI-based forecasts perform similarly to individual CI forecasts.

**Input:** Panel of LPPLS-CI time series for one bubble type/scale across all 7 or 12 countries.

**Output:** One weekly TCI series per bubble type/scale (6 series for G7; 6 for G7+BRICS).
