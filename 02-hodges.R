## Hodge's insurance premiums data ##

# Test for, and install if not found, the packages required to run this script.
pkgs <- c(
  'Matrix',
  'gamlss',
  'lme4',
  'gamlss.data',
  'TMB',
  'tidyverse',
  'mvQuad',
  'trust',
  'microbenchmark'
)
for (pkg in pkgs) {
  if (!require(pkg,character.only = TRUE,quietly = TRUE)) {
    cat(paste0("Could not find package ",pkg,", installing from CRAN.\n"))
    install.packages(pkg)
    require(pkg,character.only = TRUE,quietly = TRUE)
  }
}

data(hodges)

PLOTTEXTSIZE <- 18

globalpath <- "~/work/projects/bct-ad"
tmbpath <- file.path(globalpath,"code/02_bct_tmb")
plotpath <- file.path(globalpath,"figures")

compile(
  paste0(tmbpath,".cpp"),
  supernodal = FALSE,
  framework = "TMBad",
  libtmb = FALSE,
  tracesweep = TRUE
)
dyn.load(dynlib(tmbpath))

## Laplace approximate marginal likelihood ##
hodges_srt <- hodges[order(hodges$state), ]
lform_hodges <- lme4::lFormula(prind ~ (1|state),data=hodges_srt)
Z <- as(t(lform_hodges$reTrms$Zt),'dgTMatrix')
# Base data
tmbdata_base <- list(
  y = hodges_srt$prind,
  id = as.numeric(table(hodges_srt$state)),
  Z = Z,
  approx_type = 1L,
  re_dist = 1L,
  m = 0L,
  nn = 0,
  ww = 0
)
# Data for each candidate template
tmbdata_series4 <- tmbdata_base
tmbdata_series4$approx_type = 1
tmbdata_series4$m = 4L

tmbdata_series10 <- tmbdata_base
tmbdata_series10$approx_type = 1
tmbdata_series10$m = 10L

tmbdata_series25 <- tmbdata_base
tmbdata_series25$approx_type = 1
tmbdata_series25$m = 25L

tmbdata_series100 <- tmbdata_base
tmbdata_series100$approx_type = 1
tmbdata_series100$m = 100L

gg3 <- mvQuad::createNIGrid(1,'GLe',3)
nn3 <- mvQuad::getNodes(gg3)[ ,1]
ww3 <- mvQuad::getWeights(gg3)[ ,1]

gg9 <- mvQuad::createNIGrid(1,'GLe',9)
nn9 <- mvQuad::getNodes(gg9)[ ,1]
ww9 <- mvQuad::getWeights(gg9)[ ,1]

gg15 <- mvQuad::createNIGrid(1,'GLe',15)
nn15 <- mvQuad::getNodes(gg15)[ ,1]
ww15 <- mvQuad::getWeights(gg15)[ ,1]

gg25 <- mvQuad::createNIGrid(1,'GLe',25)
nn25 <- mvQuad::getNodes(gg25)[ ,1]
ww25 <- mvQuad::getWeights(gg25)[ ,1]

tmbdata_gl3 <- tmbdata_base
tmbdata_gl3$approx_type = 2
tmbdata_gl3$nn <- nn3
tmbdata_gl3$ww <- ww3

tmbdata_gl9 <- tmbdata_base
tmbdata_gl9$approx_type = 2
tmbdata_gl9$nn <- nn9
tmbdata_gl9$ww <- ww9

tmbdata_gl15 <- tmbdata_base
tmbdata_gl15$approx_type = 2
tmbdata_gl15$nn <- nn15
tmbdata_gl15$ww <- ww15

tmbdata_gl25 <- tmbdata_base
tmbdata_gl25$approx_type = 2
tmbdata_gl25$nn <- nn25
tmbdata_gl25$ww <- ww25

# Parameters for all templates
startingvals <- c(log(median(hodges$prind)),0,.1,0)
tmbparams <- list(
  beta = startingvals,
  theta = rep(0,4),
  u = matrix(0,nrow = length(unique(hodges_srt$state)),ncol=4)
)

# Templates
# Create 2 versions of each template:
#   1) intern = TRUE: gives the AD hessian of the Laplace approximation
#   2) intern = FALSE: makes the fitted random effects visible in the model environment
template_series4 <- MakeADFun(
  data = tmbdata_series4,
  parameters = tmbparams,
  silent = TRUE,
  DLL = "02_bct_tmb",
  random = 'u',
  intern = TRUE
)
template_series4_re <- MakeADFun(
  data = tmbdata_series4,
  parameters = tmbparams,
  silent = TRUE,
  DLL = "02_bct_tmb",
  random = 'u',
  intern = FALSE
)

template_series10 <- MakeADFun(
  data = tmbdata_series10,
  parameters = tmbparams,
  silent = TRUE,
  DLL = "02_bct_tmb",
  random = 'u',
  intern = TRUE
)
template_series10_re <- MakeADFun(
  data = tmbdata_series10,
  parameters = tmbparams,
  silent = TRUE,
  DLL = "02_bct_tmb",
  random = 'u',
  intern = FALSE
)

template_series25 <- MakeADFun(
  data = tmbdata_series25,
  parameters = tmbparams,
  silent = TRUE,
  DLL = "02_bct_tmb",
  random = 'u',
  intern = TRUE
)
template_series25_re <- MakeADFun(
  data = tmbdata_series25,
  parameters = tmbparams,
  silent = TRUE,
  DLL = "02_bct_tmb",
  random = 'u',
  intern = FALSE
)

template_series100 <- MakeADFun(
  data = tmbdata_series100,
  parameters = tmbparams,
  silent = TRUE,
  DLL = "02_bct_tmb",
  random = 'u',
  intern = TRUE
)
template_series100_re <- MakeADFun(
  data = tmbdata_series100,
  parameters = tmbparams,
  silent = TRUE,
  DLL = "02_bct_tmb",
  random = 'u',
  intern = FALSE
)


template_gl3 <- MakeADFun(
  data = tmbdata_gl3,
  parameters = tmbparams,
  silent = TRUE,
  DLL = "02_bct_tmb",
  random = 'u',
  intern = TRUE
)
template_gl3_re <- MakeADFun(
  data = tmbdata_gl3,
  parameters = tmbparams,
  silent = TRUE,
  DLL = "02_bct_tmb",
  random = 'u',
  intern = FALSE
)

template_gl9 <- MakeADFun(
  data = tmbdata_gl9,
  parameters = tmbparams,
  silent = TRUE,
  DLL = "02_bct_tmb",
  random = 'u',
  intern = TRUE
)
template_gl9_re <- MakeADFun(
  data = tmbdata_gl9,
  parameters = tmbparams,
  silent = TRUE,
  DLL = "02_bct_tmb",
  random = 'u',
  intern = FALSE
)

template_gl15 <- MakeADFun(
  data = tmbdata_gl15,
  parameters = tmbparams,
  silent = TRUE,
  DLL = "02_bct_tmb",
  random = 'u',
  intern = TRUE
)
template_gl15_re <- MakeADFun(
  data = tmbdata_gl15,
  parameters = tmbparams,
  silent = TRUE,
  DLL = "02_bct_tmb",
  random = 'u',
  intern = FALSE
)

template_gl25 <- MakeADFun(
  data = tmbdata_gl25,
  parameters = tmbparams,
  silent = TRUE,
  DLL = "02_bct_tmb",
  random = 'u',
  intern = TRUE
)
template_gl25_re <- MakeADFun(
  data = tmbdata_gl25,
  parameters = tmbparams,
  silent = TRUE,
  DLL = "02_bct_tmb",
  random = 'u',
  intern = FALSE
)
## Optimization ##

fit_gamlssmod <- function() {
  gamlss(
    prind ~ random(factor(state)),
    sigma.formula = ~random(factor(state)),
    nu.formula = ~random(factor(state)),
    tau.formula = ~random(factor(state)),
    family = BCT,
    data = hodges
  )
}
tm <- Sys.time()
gamlssmod <- fit_gamlssmod()
time_gamlss <- as.numeric(difftime(Sys.time(),tm,units = 'secs'))

optimize_model_trd <- function(template) {
  objfun <- function(p) with(template,list(value = fn(p),gradient = as.numeric(gr(p)),hessian = he(p)))
  trust::trust(objfun,template$par,2,5)
}

# Series
tm <- Sys.time()
opt_series4 <- optimize_model_trd(template_series4)
time_series4 <- as.numeric(difftime(Sys.time(),tm,units = 'secs'))

tm <- Sys.time()
opt_series10 <- optimize_model_trd(template_series10)
time_series10 <- as.numeric(difftime(Sys.time(),tm,units = 'secs'))

tm <- Sys.time()
opt_series25 <- optimize_model_trd(template_series25)
time_series25 <- as.numeric(difftime(Sys.time(),tm,units = 'secs'))

tm <- Sys.time()
opt_series100 <- optimize_model_trd(template_series100)
time_series100 <- as.numeric(difftime(Sys.time(),tm,units = 'secs'))

# GLQ
tm <- Sys.time()
opt_gl3 <- optimize_model_trd(template_gl3)
time_gl3 <- as.numeric(difftime(Sys.time(),tm,units = 'secs'))

tm <- Sys.time()
opt_gl9 <- optimize_model_trd(template_gl9)
time_gl9 <- as.numeric(difftime(Sys.time(),tm,units = 'secs'))

tm <- Sys.time()
opt_gl15 <- optimize_model_trd(template_gl15)
time_gl15 <- as.numeric(difftime(Sys.time(),tm,units = 'secs'))

tm <- Sys.time()
opt_gl25 <- optimize_model_trd(template_gl25)
time_gl25 <- as.numeric(difftime(Sys.time(),tm,units = 'secs'))

# More thorough timing
opt_times <- microbenchmark::microbenchmark(
  fit_gamlssmod(),

  optimize_model_trd(template_series4),
  optimize_model_trd(template_series10),
  optimize_model_trd(template_series25),
  optimize_model_trd(template_series100),

  optimize_model_trd(template_gl3),
  optimize_model_trd(template_gl9),
  optimize_model_trd(template_gl15),
  optimize_model_trd(template_gl25),
  times = 100
)
opt_times_summary <- as.data.frame(opt_times) %>%
  group_by(expr) %>%
  summarize(median = median(time))
gamlss_median <- opt_times_summary[opt_times_summary$expr=='fit_gamlssmod()','median',drop=TRUE]
opt_times_summary_ratios <- opt_times_summary[opt_times_summary$expr!='fit_gamlssmod()', ]
opt_times_summary_ratios$ratio <- 100 * opt_times_summary_ratios$median / gamlss_median

write_csv(opt_times_summary_ratios,file = file.path(plotpath,"hodges_opttimes.csv"))

# Compare solutions
opt_coefs <- rbind(
  opt_series4$argument,
  opt_series10$argument,
  opt_series25$argument,
  opt_series100$argument,
  
  opt_gl3$argument,
  opt_gl9$argument,
  opt_gl15$argument,
  opt_gl25$argument
)
# get the sd from the log sd
opt_coefs[ ,5:8] <- exp(opt_coefs[ ,5:8])
colnames(opt_coefs) <- c("beta_mu","beta_sigma","beta_nu","beta_tau",paste0("alpha",1:4))
opt_coefs <- as_tibble(opt_coefs) %>%
  add_column(method = c(rep("series",4),rep("gl",4)),tuning = as.character(c(4,10,25,100,3,9,15,25)))
# add the one from gamlss
gamlss_betas <- as.numeric(summary(gamlssmod)[ ,1])
gamlss_alphas <- c(
  getSmo(gamlssmod,what = 'mu')$sigb,
  getSmo(gamlssmod,what = 'sigma')$sigb,
  getSmo(gamlssmod,what = 'nu')$sigb,
  getSmo(gamlssmod,what = 'tau')$sigb
)
opt_coefs <- opt_coefs %>% add_row(
  beta_mu = gamlss_betas[1],
  beta_sigma = gamlss_betas[2],
  beta_nu = gamlss_betas[3],
  beta_tau = gamlss_betas[4],
  
  alpha1 = gamlss_alphas[1],
  alpha2 = gamlss_alphas[2],
  alpha3 = gamlss_alphas[3],
  alpha4 = gamlss_alphas[4],
  
  method = "gamlss",
  tuning = "0"
)
write_csv(opt_coefs,file = file.path(plotpath,"hodges_optcoef.csv"))
# note: results for mu aren't comparable to gamlss since this was done on the natural scale in gamlss, but the log scale in TMB
