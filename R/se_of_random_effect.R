library(merDeriv)
library(lme4)
library(broom.mixed)

data = gen_two_level_data(m = 100, n = 5, sigma_b2 = 2.5, seed = 1083)

outp <- gen_mor_estimate(m = 100, n = 50, sigma_b2 = 2.5, seed = 1083)

vv <- merDeriv::vcov.glmerMod(outp$model, full = TRUE, ranpar = "sd")

sqrt(diag(vv)[4])



# ===============

m1 <- lmer(Reaction ~ Days + (Days|Subject), data = sleepstudy)

vv1 <- vcov(m1, full=TRUE)

sqrt(diag(vv1))
