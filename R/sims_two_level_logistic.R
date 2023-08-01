m = 50 # number of cluster
n = 5 # size of each cluster
N = m*n 

set.seed(1083)

x1c <- rnorm(N)
x2 <- rnorm(N)
x2b <- ifelse(x2 <= 0.5, 0, 1)

sigma_b2 <- 2.5
beta0 <- rnorm(1)
beta1 <- rnorm(1)
beta2 <- rnorm(1)

uj <- rep(rnorm(m, 0, sqrt(sigma_b2)), each = n)
eta_ij <- beta0 + beta1*x1c + beta2*x2b + uj
pi_ij <- exp(eta_ij) / (1 + exp(eta_ij))

yij <- rbinom(N, size = 1, prob = pi_ij)

multi_data <- data.frame(
  cluster = rep(seq(m), each = n),
  id = rep(seq(n), times = m),
  uj = uj,
  X1c = x1c,
  X2b= x2b,
  Yij = yij
)

library(lme4)
library(broom.mixed)

multi_model <- glmer(Yij ~ X1c + X2b + (1 | cluster), 
                     family = binomial("logit"),
                     data = multi_data)

summary(multi_model)
output_df <- tidy(multi_model)

sigma_b2_hat <- as.numeric(output_df[output_df$effect == "ran_pars", "estimate"]) ^ 2

MOR <- exp(sqrt(2 * sigma_b2) * qnorm(0.75))
MOR_hat <- exp(sqrt(2 * sigma_b2_hat) * qnorm(0.75))
