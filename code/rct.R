library(DeclareDesign)

declaration_18.1 <-
  declare_model(N = 100,
                U = rnorm(N),
                potential_outcomes(Y ~ 0.2 * Z + U)) +
  declare_inquiry(ATE = mean(Y_Z_1 - Y_Z_0)) +
  declare_assignment(Z = complete_ra(N, prob = 0.5)) +
  declare_measurement(Y = reveal_outcomes(Y ~ Z)) +
  declare_estimator(Y ~ Z, inquiry = "ATE")

diagnosis_18.1 <- diagnose_design(declaration_18.1)

library(marginaleffects)

tidy_slope <- function(x) {
  tidy(avg_slopes(x, data = x$data))
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


diagnosis_11.5 <- diagnose_design(declaration_11.5)


N <- 1000
age <- runif(min = 22, max = 28, n = N)
quality <- rbinom(prob = 0.5, size = N)
married <- rbinom(n = 1, prob = 0.5, size = N)
woman <- rbinom(n = 1, prob = 0.5, size = N)

sim_d <- data.frame(age = age,
  quality = quality, married = married, woman = woman)

