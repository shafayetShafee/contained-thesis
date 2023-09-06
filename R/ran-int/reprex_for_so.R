# copied from
# http://rstudio-pubs-static.s3.amazonaws.com/28864_dd1f084207d54f5ea67c9d1a9c845d01.html
vcov.VarCorr.merMod <- function(object,fit,...) {
  if (isREML(fit)) {
    warning("refitting model with ML")
    fit <- refitML(fit)
  }
  if (!require("numDeriv")) stop("numDeriv package required")
  useSc <- attr(object,"useSc")
  dd <- lme4:::devfun2(fit,useSc=useSc,signames=FALSE)
  vdd <- as.data.frame(object,order="lower.tri")
  pars <- vdd[,"sdcor"]
  npar0 <- length(pars)
  if (isGLMM(fit)) {
    pars <- c(pars,fixef(fit))
  }
  hh1 <- hessian(dd,pars)
  vv2 <- 2*solve(hh1)
  if (isGLMM(fit)) {
    vv2 <- vv2[1:npar0,1:npar0,drop=FALSE]
  }
  nms <- apply(vdd[,1:3],1,
               function(x) paste(na.omit(x),collapse="."))
  dimnames(vv2) <- list(nms,nms)
  return(vv2)
}

library(lme4)

l = 20 # for level 3
m = 10 # for level 2
n = 15 # for level 1
N = l*m*n 

set.seed(344920)

x1c <- rnorm(N)
x2 <- rnorm(N)
x2b <- ifelse(x2 <= 0.5, 0, 1)

sigma_u_sq <- 2
sigma_v_sq <- 2.5

beta0 <- -1.85
beta1 <- 1.75
beta2 <- 0.67

ujk <- rep(rnorm(l*m, 0, sqrt(sigma_u_sq)), each = n)
vk <- rep(rnorm(l, 0, sqrt(sigma_v_sq)), each = m*n)

eta_ijk <- beta0 + beta1*x1c + beta2*x2b + ujk + vk
pi_ijk <- exp(eta_ijk) / (1 + exp(eta_ijk))
yijk <- rbinom(N, size = 1, prob = pi_ijk)

multi_data_int <- data.frame(
  ea = rep(seq(l), each = m*n),
  hh = rep(rep(seq(m), each = n), times = l),
  id = rep(seq(n), times = l*m),
  ujk = ujk,
  vk = vk,
  X1c = x1c,
  X2b = x2b,
  Yijk = yijk
)

head(multi_data_int, 20)

model <- glmer(Yijk ~ X1c + X2b + (1 | ea) + (1 | ea:hh),
               data = multi_data_int, family = "binomial")

s = summary(model)

vcov(VarCorr(model), model)

