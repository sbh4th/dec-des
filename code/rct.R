library(tidyverse)
library(here)
library(DeclareDesign)
library(marginaleffects)
library(tinytable)
library(modelsummary)

declaration_18.1 <-
  declare_model(N = 100,
                U = rnorm(N),
                potential_outcomes(Y ~ 0.2 * Z + U)) +
  declare_inquiry(ATE = mean(Y_Z_1 - Y_Z_0)) +
  declare_assignment(Z = complete_ra(N, prob = 0.5)) +
  declare_measurement(Y = reveal_outcomes(Y ~ Z)) +
  declare_estimator(Y ~ Z, inquiry = "ATE") +
  declare_estimator(Z ~ U, label = "Balance")

dd_data <- run_design(declaration_18.1)

diagnosis_18.1 <- diagnose_design(declaration_18.1)

library(marginaleffects)

tidy_slope <- function(x) {
  tidy(avg_slopes(x, var = x$Z))
}

N <- 10

declaration_11.5 <-
  declare_model(N = N,
                U = rnorm(N),
                potential_outcomes(Y ~ rbinom(N, 1, prob = 0.2 * Z + 0.6))) +
  declare_inquiry(ATE = mean(Y_Z_1 - Y_Z_0)) +
  declare_assignment(Z = complete_ra(N, prob = 0.5)) +
  declare_measurement(Y = reveal_outcomes(Y ~ Z)) +
  declare_estimator(Y ~ Z,
                    inquiry = "ATE",
                    term = "Z",
                    label = "OLS") +
  declare_estimator(
    Y ~ Z,
    .method = glm,
    family = binomial("logit"),
    .summary = tidy_slope,
    inquiry = "ATE",
    term = "Z",
    label = "logit"
  ) +
  declare_estimator(
    Y ~ Z,
    .method = glm,
    family = binomial("probit"),
    .summary = tidy_slope,
    inquiry = "ATE",
    term = "Z",
    label = "probit"
  ) 


diagnosis_11.5 <- diagnose_design(declaration_11.5, sims = 30)


N <- 1e6
age <- runif(min = 22, max = 28, n = N)
quality <- rbinom(N, 1, prob = 0.5)
married <- rbinom(N, 1, prob = 0.5)
woman <- rbinom(N, 1, prob = 0.5)

prob_to_logit_coef <- function(p0, delta) {
  # p0    = baseline probability (quality=0)
  # delta = desired change in probability when covariate=1
  qlogis(p0 + delta) - qlogis(p0)
}

# Baseline probability when all covariates = 0
p0 <- 0.10

# Desired probability changes:
delta_quality <- 0.20   # add 0.25 when quality=1
delta_married <- 0.10   # add 0.10 when married=1
delta_woman   <- 0.05   # add 0.05 when woman=1

# Convert to logit coefficients
alpha <- qlogis(p0)
beta_quality <- prob_to_logit_coef(p0, delta_quality)
beta_married <- prob_to_logit_coef(p0, delta_married)
beta_woman   <- prob_to_logit_coef(p0, delta_woman)

# Simulate
sim_d <- data.frame(age = age,
                    quality = quality,
                    married = married,
                    woman = woman) %>%
  mutate(
    linpred = alpha + (beta_quality * quality) +
      (beta_married * married) + (beta_woman * woman),
    p_lpm = pmin(pmax(linpred, 0), 1),
    linprob = alpha + (delta_quality * quality) +
      (delta_married * married) + (delta_woman * woman),
    p = plogis(linpred),
    cb = rbinom(n(), 1, prob = p),
    y = rbinom(n(), 1, prob = p_lpm))


modelsummary::datasummary_skim(sim_d)

mod <- glm(
  cb ~ quality + married + woman,
  data = sim_d, family = "binomial")

mod_lmp <- lm_robust(
  cb ~ quality + married + woman,
  data = sim_d)

baseline <- avg_predictions(mod,
  newdata = datagrid(quality = 0,
  married = 0, woman = 0))
tidy(baseline)

me <- avg_slopes(mod)
tidy(me)



N <- 1e6
# desired probability-scale parameters
p0 <- 0.10              # baseline P(Y=1) when all covariates = 0
delta_quality <- 0.5   # additive probability effect of quality (treated = 1)
delta_married <- 0.4
delta_woman <- 0.3

# covariate generation (tweak probs as desired)
quality  <- rbinom(N, 1, 0.5)
married  <- rbinom(N, 1, 0.5)
woman    <- rbinom(N, 1, 0.5)

# build the conditional probability directly on probability scale:
p <- p0 + delta_quality * quality + delta_married * married + delta_woman * woman

# OPTIONAL: keep probabilities strictly in [0,1]
p <- pmin(pmax(p, 0), 1)

# draw binary outcome
y <- rbinom(N, 1, prob = p)

# data.frame and OLS
dat <- data.frame(y, quality, married, woman)

ols <- lm(y ~ quality + married + woman, data = dat)
lmr <- lm_robust(y ~ quality + married + woman, data = dat)
modelsummary(list("OLS" = ols, "Robust" = lmr))


M <- 
  declare_model(
    N = 100, 
    potential_outcomes(Y ~ rbinom(N, size = 1, prob = 0.1 * Z + 0.5)),
    Z = ifelse()
  )

M()

