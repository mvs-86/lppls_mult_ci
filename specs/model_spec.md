# Model Specification

Mathematical model specifications extracted from:
> Gabauer, Gupta, Karmakar & Nielsen (2022). *Stock Market Bubbles and the Forecastability of Gold Returns (and Volatility)*. University of Pretoria Working Paper 2022-28.

All equation numbers refer to the paper.

---

## 1. LPPLS Model (§2.1)

### 1.1 Original formulation (Eq. 1)

$$E[\ln p(t)] = A + B(t_c - t)^m + C(t_c - t)^m \cos\!\bigl(\omega \ln(t_c - t)^m - \phi\bigr)$$

**Variable definitions:**

| Symbol | Definition |
|--------|-----------|
| $p(t)$ | Log price-dividend ratio at time $t$ |
| $t_c$ | Critical time — date of bubble termination |
| $A$ | Expected log-value of the series at $t_c$ |
| $B$ | Amplitude of power-law acceleration ($B < 0$ for positive bubble, $B > 0$ for negative bubble) |
| $C$ | Relative magnitude of log-periodic oscillations |
| $m$ | Exponent of power-law growth |
| $\omega$ | Frequency of log-periodic oscillations |
| $\phi$ | Phase shift parameter |

### 1.2 Filimonov–Sornette reparameterisation (Eq. 2)

Phase $\phi$ is eliminated and $C$ is expanded as $C_1 = C\cos\phi$, $C_2 = C\cos\phi$:

$$E[\ln p(t)] = A + B \cdot f + C_1 \cdot g + C_2 \cdot h$$

where:

$$f = (t_c - t)^m$$
$$g = (t_c - t)^m \cos[\omega \ln(t_c - t)]$$
$$h = (t_c - t)^m \sin[\omega \ln(t_c - t)]$$

**Parameters:** 3 nonlinear $\{t_c, m, \omega\}$ and 4 linear $\{A, B, C_1, C_2\}$.

### 1.3 Objective function — sum of squared residuals (Eq. 3)

$$F(t_c, m, \omega, A, B, C_1, C_2) = \sum_{i=1}^{N} \bigl[\ln p(\tau_i) - A - B(f_i) - C_1(g_i) - C_2(h_i)\bigr]^2$$

### 1.4 Concentrated cost function (Eq. 4)

$$F_1(t_c, m, \omega) = \min_{A,B,C_1,C_2} F(t_c, m, \omega, A, B, C_1, C_2) = F(t_c, m, \omega, \hat{A}, \hat{B}, \hat{C}_1, \hat{C}_2)$$

### 1.5 Linear parameters — analytical solution (Eq. 5 & 6)

$$\{\hat{A}, \hat{B}, \hat{C}_1, \hat{C}_2\} = \arg\min_{A,B,C_1,C_2} F(t_c, m, \omega, A, B, C_1, C_2)$$

Solved via the normal equations:

$$\begin{pmatrix} N & \sum f_i & \sum g_i & \sum h_i \\ \sum f_i & \sum f_i^2 & \sum f_i g_i & \sum f_i h_i \\ \sum g_i & \sum f_i g_i & \sum g_i^2 & \sum g_i h_i \\ \sum h_i & \sum f_i h_i & \sum g_i h_i & \sum h_i^2 \end{pmatrix} \begin{pmatrix} \hat{A} \\ \hat{B} \\ \hat{C}_1 \\ \hat{C}_2 \end{pmatrix} = \begin{pmatrix} \sum \ln p_i \\ \sum f_i \ln p_i \\ \sum g_i \ln p_i \\ \sum h_i \ln p_i \end{pmatrix}$$

### 1.6 Nonlinear optimisation (Eq. 7)

$$\{\hat{t}_c, \hat{m}, \hat{\omega}\} = \arg\min_{t_c,m,\omega} F_1(t_c, m, \omega)$$

Solved via **SLSQP** (Sequential Least Squares Programming; Kraft 1988). In R: `nloptr()` with algorithm `"NLOPT_LD_SLSQP"`.

---

## 2. Multi-Scale LPPLS Confidence Indicator (§2.2)

### 2.1 Window construction

For each weekly observation $t_2$, sweep nested windows $[t_1, t_2]$:

- Time series sampled in steps of **5 trading days**.
- Windows iterated in steps of **2 trading days**.

| Scale | Window length $dt = t_2 - t_1$ | Number of fits |
|-------|-------------------------------|----------------|
| Short | $[30, 90]$ trading days | $(90 - 30)/2 = 30$ |
| Medium | $[90, 300]$ trading days | $(300 - 90)/2 = 105$ |
| Long | $[300, 745]$ trading days | $(745 - 300)/2 = 223$ |

### 2.2 Bubble classification

- **Positive bubble** (crash precursor): qualified fit has $\hat{B} < 0$
- **Negative bubble** (rally precursor): qualified fit has $\hat{B} > 0$

### 2.3 Filter conditions

A fit is **qualified** only if all of the following hold simultaneously:

$$m \in [0.01,\ 0.99]$$

$$\omega \in [2,\ 15]$$

$$t_c \in \bigl[\max(t_2 - 60,\ t_2 - 0.5 \cdot dt),\ \min(252,\ t_2 + 0.5 \cdot dt)\bigr]$$

$$O > 2.5 \qquad \text{where } O = \frac{\omega}{2\pi} \ln\!\left(\frac{t_c - t_1}{t_c - t_2}\right)$$

$$D > 0.5 \qquad \text{where } D = \frac{m|B|}{\omega|C|},\quad |C| = \sqrt{C_1^2 + C_2^2}$$

**$O$** measures the number of log-periodic oscillations in the window.
**$D$** measures the relative damping of oscillations to the trend.

### 2.4 Confidence indicator

$$\text{LPPLS-CI}_{\text{scale}} = \frac{\text{number of qualified fits at scale}}{\text{total fits at scale}} \in [0, 1]$$

This yields **6 indicators per country per week**: positive and negative × short, medium, long.

---

## 3. Forecasting Model (§2.3)

### 3.1 Model specifications

**Model 1** — AR(1) benchmark:

$$y_t = a_0 + a_1 y_{t-1} + e_t$$

**Model 2** — AR(1) augmented with LPPLS-CI predictors:

$$y_t = a_0 + a_1 y_{t-1} + \sum_{j=0}^{p} \beta_j x_{t,j} + e_t$$

where $y_t$ is the weekly gold log-return (or TVPGARCH volatility), $x_{t,j}$ are the LPPLS-CI predictors ($p = 42$ for G7 only; $p = 72$ for G7 + BRICS).

### 3.2 LASSO estimation of Model 2 (Step 1)

$$\hat{\kappa} = \arg\min_{\kappa} \left( \sum_{t=1}^{T} \left(y_t - a_0 - a_1 y_{t-1} - \sum_{j=1}^{p} \beta_j x_{t,j}\right)^2 + \lambda\left(|a_0| + |a_1| + \sum_j |\beta_j|\right) \right)$$

where $\kappa = (a_0, a_1, \beta_1, \ldots, \beta_p)$.

- $\lambda$ selected by cross-validation via `cv.glmnet(..., lambda = "lambda.min")`.
- $a_0$ and $a_1$ are also shrunk (included in the penalty).

### 3.3 Block-mean residual forecast procedure

**Step 1.** Obtain LASSO residuals:

$$\hat{e}_i = y_i - \hat{y}_i$$

**Step 2.** Form block means of residuals over horizon $h$:

$$\tilde{e}_{i,h} = \frac{\hat{e}_1 + \cdots + \hat{e}_h}{h}$$

**Step 3.** Point forecast of $(y_{T+1} + \cdots + y_{T+h})/h$:

$$\hat{F}_h = \frac{\hat{y}_{T+1} + \cdots + \hat{y}_{T+h}}{h} + \operatorname{mean}(\tilde{e}_h)$$

### 3.4 Forecast evaluation

**Clark–West (CW) test** (Clark & West 2007) — for nested model comparison (Model 2 vs. Model 1).

- Pseudo-out-of-sample, rolling window with training lengths $T \in \{100, 250, 500, 1000\}$.
- Forecast horizons $h \in \{1, 2, 4, 8, 12\}$ weeks.

**Diebold–Mariano (DM) test** (Diebold & Mariano 1995) — for non-nested comparison (positive-only vs. negative-only LPPLS-CI models).

---

## 4. TVPGARCH Volatility Model (§2.4)

### 4.1 Model specification (Eq. 8)

$$y_t \sim N(0, \sigma_t^2) \quad \text{with} \quad \sigma_t^2 = \alpha_0(t/n) + \alpha_1(t/n)\,y_{t-1}^2 + \beta_1(t/n)\,\sigma_{t-1}^2$$

**Variable definitions:**

| Symbol | Definition |
|--------|-----------|
| $y_t$ | Weekly gold log-return |
| $\sigma_t^2$ | Conditional variance at time $t$ |
| $\alpha_0(\cdot)$ | Time-varying intercept function |
| $\alpha_1(\cdot)$ | Time-varying ARCH coefficient function |
| $\beta_1(\cdot)$ | Time-varying GARCH coefficient function |
| $\theta = (\alpha_0, \alpha_1, \beta_1)$ | Parameter vector |

### 4.2 Kernel-based local likelihood estimation (Eq. 9)

$$\hat{\theta}_{b_n}(t) = \operatorname{argmin}_{\theta \in \Theta} \sum_{i=1}^{n} K\!\left(\frac{t - i/n}{b_n}\right) \ell(y_i, X_i, \theta), \quad t \in [0, 1]$$

- $K(\cdot)$: **Epanechnikov kernel**, $K(u) = \tfrac{3}{4}(1 - u^2)\mathbf{1}(|u| \le 1)$
- $b_n \in [0, 1]$: bandwidth (default $b_n = 0.25$; other values available upon request)
- $X_i = y_{i-1}$ (univariate GARCH)

### 4.3 Quasi log-likelihood (Eq. 9, expanded)

$$\ell(y_i, X_i, \theta') = -\frac{1}{2}\log(\sigma_i^2) + \frac{y_i^2}{\sigma_i^2} \quad \text{with} \quad \sigma_i^2 = \alpha_0 + \alpha_1 y_{i-1}^2 + \beta_1 \sigma_{i-1}^2$$

The estimated $\hat{\sigma}^2$ series replaces $y_t$ as the dependent variable in the forecasting model (§3) to produce volatility forecasts.

---

## 5. Total Connectedness Index via TVP-VAR (Appendix)

Used as an alternative aggregation of the 6 per-country LPPLS-CI indicators into a single connectedness measure (TCI). Applied separately to the G7 (6 bubble indicators) and G7+BRICS (6 bubble indicators).

### 5.1 TVP-VAR system (Eqs. A.1–A.2)

$$\mathbf{z}_t = \mathbf{B}_t \mathbf{z}_{t-1} + \mathbf{u}_t, \qquad \mathbf{u}_t \sim N(\mathbf{0}, \mathbf{S}_t)$$

$$\operatorname{vec}(\mathbf{B}_t) = \operatorname{vec}(\mathbf{B}_{t-1}) + \mathbf{v}_t, \qquad \mathbf{v}_t \sim N(\mathbf{0}, \mathbf{R}_t)$$

**Variable definitions:**

| Symbol | Dimension | Definition |
|--------|-----------|-----------|
| $\mathbf{z}_t$ | $k \times 1$ | Vector of bubble indicators at time $t$ |
| $\mathbf{B}_t$ | $k \times k$ | Time-varying VAR coefficient matrix |
| $\mathbf{S}_t$ | $k \times k$ | Time-varying variance-covariance matrix |
| $\mathbf{v}_t$ | $k^2 \times 1$ | State innovation vector |
| $\mathbf{R}_t$ | $k^2 \times k^2$ | State innovation covariance matrix |

Lag length $p = 1$ (selected by BIC).

### 5.2 Generalised Forecast Error Variance Decomposition (GFEVD)

TVP-VAR is transformed to its TVP-VMA representation:

$$\mathbf{z}_t = \sum_{s=1}^{p} \mathbf{B}_t \mathbf{z}_{t-s} + \mathbf{u}_t = \sum_{j=0}^{\infty} \mathbf{A}_{jt} \mathbf{u}_{t-j}$$

The GFEVD (unnormalised and normalised):

$$\psi^g_{ij,t}(H) = \frac{S^{-1}_{ii,t} \sum_{t=1}^{H-1} (\boldsymbol{\iota}_i' \mathbf{A}_t \mathbf{S}_t \boldsymbol{\iota}_j)^2}{\sum_{k=1}^{K} \sum_{t=1}^{H-1} (\boldsymbol{\iota}_i \mathbf{A}_t \mathbf{S}_t \mathbf{A}_t' \boldsymbol{\iota}_i)}, \qquad \tilde{\psi}^g_{ij,t}(H) = \frac{\psi^g_{ij,t}(H)}{\sum_{k=1}^{K} \phi^g_{ij,t}(H)}$$

where $\boldsymbol{\iota}_i$ is a zero vector with unity at position $i$, and $\sum_{j=1}^k \tilde{\psi}^g_{ij,t}(H) = 1$, $\sum_{j=1}^k \tilde{\psi}^g_{ij,t}(H) = k$.

### 5.3 Total Connectedness Index (Eq. A.3)

$$TCI = \frac{\sum_{i,j=1,\, i\neq j}^{k} \tilde{\psi}^g_{ij,t}(H)}{\sum_{i,j=1}^{k} \tilde{\psi}^g_{ij,t}(H)}, \qquad 0 \le C^g_t(H) \le 1$$

The TCI series (one per scale: short/medium/long, separately for G7 and G7+BRICS) is then used as a compact alternative to the full set of individual LPPLS-CI predictors in the forecasting model.
