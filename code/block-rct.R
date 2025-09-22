library(tidyverse)
library(here)
library(DeclareDesign)
library(scales)

declaration_18.4 <-
  declare_model(
    N = 500,
    X = rep(c(0, 1), each = N / 2),
    U = rnorm(N, sd = 0.25),
    potential_outcomes(Y ~ 0.2 * Z + X + U)
  ) +
  declare_assignment(
    Z = block_ra(blocks = X, block_prob = c(0.2, 0.5)),
    probs =
      obtain_condition_probabilities(
        assignment = Z, 
        blocks = X, 
        block_prob = c(0.2, 0.5)),
    ipw = 1 / probs) +
  declare_measurement(Y = reveal_outcomes(Y ~ Z)) +
  declare_estimator(
    Y ~ Z,
    covariates = ~ X,
    .method = lm_lin,
    weights = ipw,
    label = "Lin"
  ) +
  declare_estimator(
    Y ~ X + Z,
    .method = lm_robust,
    weights = ipw,
    label = "Robust"
  )

diagnosis_18.4 <- diagnose_design(declaration_18.4)
diagnosis_18.4


set.seed(343)
fixed_pop <-
  fabricate(
    N = 12,
    X = rep(c(0, 1), each = N/2),
    U = rnorm(N, sd = 0.25),
    potential_outcomes(Y ~ 0.2*Z + X + U)
  )

ra_declaration <- declare_ra(N = 12)
permutation_matrix <- obtain_permutation_matrix(ra_declaration)
estimates_list <- list()


gg_df <-
  diagnosis_18.4 |>
  dplyr::select(sim_ID, estimator, estimate) |>
  pivot_wider(names_from = estimator, values_from = estimate) |>
  mutate(balanced = balance == 0)


gg_df <-
  diagnosis_18.4 |>
  select(sim_ID, estimator, estimate) |>
  pivot_wider(names_from = estimator, values_from = estimate) |>
  mutate(balanced = balance == 0)

g <-
  ggplot(gg_df) +
  aes(DIM, fill = balanced) +
  geom_histogram(
    aes(y = ..count.. / sum(..count..)),
    position = "identity",
    binwidth = 0.11,
    alpha = 0.7
  ) +
  annotate(
    "text",
    x = 0.4,
    y = .04,
    label = "Assignments that\nexactly balance X",
    color = dd_palette("dd_dark_blue"),
    hjust = 0
  ) +
  annotate(
    "text",
    x = -0.1,
    y = .06,
    label = "Assignments that do not\nexactly balance X",
    color = dd_palette("dd_gray"),
    hjust = 1
  ) +
  scale_fill_manual(values = rev(dd_palette("two_color_gray"))) +
  # scale_color_manual(values = dd_palette("two_color_gray")) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     breaks = seq(0, 0.1, 0.02)) +
  theme_dd() +
  theme(axis.title = element_text(size = 9)) +
  labs(x = "Simulated effect estimate",
       y = "Percent of possible random assignments")

g



declaration <- declare_ra(N=100, m_each=c(50, 50))
declaration

Z <- conduct_ra(declaration)
table(Z)

probs <- obtain_condition_probabilities(declaration, Z)
table(probs, Z)


# Blocked and Clustered Random Assignment Declarations

clusters <- rep(letters, times=1:26)
blocks <- rep(NA, length(clusters))
blocks[clusters %in% letters[1:5]] <- "block_1"
blocks[clusters %in% letters[6:10]] <- "block_2"
blocks[clusters %in% letters[11:15]] <- "block_3"
blocks[clusters %in% letters[16:20]] <- "block_4"
blocks[clusters %in% letters[21:26]] <- "block_5"

table(blocks, clusters)

declare_ra(clusters = clusters, blocks = blocks)
declare_ra(clusters = clusters, blocks = blocks, prob_each = c(.2, .5, .3))


library(randomizr)
# number of jobs
N <- 1000
age <- round(runif(min = 22, max = 28, n = N),0)
quality <- rbinom(N, 1, prob = 0.5)
married <- rbinom(N, 1, prob = 0.5)
woman <- rbinom(N, 1, prob = 0.5)

# create dataset
sim_d <- data.frame(job_id = 1:N,
  age = age, quality = quality,
  married = married, woman = woman)

d <- data.frame(job_id = 1:N)


bernoulli_logit <- fabricate(
  N = N, x = 10 * rnorm(N),
  binary = draw_binary(latent = x, link = "probit")
)
