source(here::here("R/sim_funcs.R"))

# problematic seed m = 10, n = 5
# 869455 652503
gen_mor_estimate(m = 10, n = 5, sigma_b2 = 2.5, seed = 869455)
gen_mor_estimate(m = 10, n = 5, sigma_b2 = 2.5, seed = 652503)


# problematic seed m = 10, n = 10
# 980506
gen_mor_estimate(m = 10, n = 10, sigma_b2 = 2.5, seed = 980506)
