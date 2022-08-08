# Alternative models for Hodge's insurance premium data #


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
PLOTTEXTSIZESMALL <- 12

## SET THESE PATHS ##

globalpath <- "~/work/projects/bct-ad"
tmbpath <- file.path(globalpath,"code/03_bct_gaussian")
tmbpath2 <- file.path(globalpath,"code/03_bct_heavytail")
plotpath <- file.path(globalpath,"figures")

## END SET THESE PATHS ##

compile(
  paste0(tmbpath,".cpp"),
  supernodal = FALSE,
  framework = "TMBad",
  libtmb = FALSE,
  tracesweep = TRUE
)
compile(
  paste0(tmbpath2,".cpp"),
  supernodal = FALSE,
  framework = "TMBad",
  libtmb = FALSE,
  tracesweep = TRUE
)
dyn.load(dynlib(tmbpath))
dyn.load(dynlib(tmbpath2))

## Laplace approximate marginal likelihood ##
hodges_srt <- hodges[order(hodges$state), ]
lform_hodges <- lme4::lFormula(prind ~ (1|state),data=hodges_srt)
Z <- as(t(lform_hodges$reTrms$Zt),'dgTMatrix')

gg15 <- mvQuad::createNIGrid(1,'GLe',15)
nn15 <- mvQuad::getNodes(gg15)[ ,1]
ww15 <- mvQuad::getWeights(gg15)[ ,1]

tmbdata <- list(
  y = hodges_srt$prind,
  id = as.numeric(table(hodges_srt$state)),
  Z = Z,
  approx_type = 2L,
  re_dist = 1L,
  m = 0L,
  nn = nn15,
  ww = ww15
)
startingvals <- c(log(median(hodges$prind)),0,.1,0)
tmbparams <- list(
  beta = startingvals,
  theta = rep(0,2),
  u = matrix(0,nrow = length(unique(hodges_srt$state)),ncol=2)
)

template_gaussian <- MakeADFun(
  data = tmbdata,
  parameters = tmbparams,
  silent = TRUE,
  DLL = "03_bct_gaussian",
  random = 'u',
  intern = TRUE
)
template_gaussian_re <- MakeADFun(
  data = tmbdata,
  parameters = tmbparams,
  silent = TRUE,
  DLL = "03_bct_gaussian",
  random = 'u',
  intern = FALSE
)

# Fit the Gaussian model
objfun <- function(p) with(template_gaussian,list(value = fn(p),gradient = as.numeric(gr(p)),hessian = he(p)))
opt_gaussian <- trust::trust(objfun,template_gaussian$par,2,5)
# coefficients
betas_gaussian <- opt_gaussian$argument[1:4]
alphas_gaussian <- exp(opt_gaussian$argument[5:6])

# random effects
template_gaussian_re$fn(opt_gaussian$argument)
re_gaussian <- matrix(with(template_gaussian_re$env,last.par.best[random]),ncol=2)

# predicted values
# point predictions- medians
state_point_predictions <- exp(betas_gaussian[1] + re_gaussian[ ,1])
state_observations <- with(hodges_srt,tapply(prind,state,median))

## heavy-tailed random effects ##

tmbparams_heavy <- tmbparams
tmbdata_heavy <- tmbdata

template_heavytail <- MakeADFun(
  data = tmbdata_heavy,
  parameters = tmbparams_heavy,
  silent = TRUE,
  DLL = "03_bct_heavytail",
  random = 'u',
  intern = TRUE
)
template_heavytail_re <- MakeADFun(
  data = tmbdata_heavy,
  parameters = tmbparams_heavy,
  silent = TRUE,
  DLL = "03_bct_heavytail",
  random = 'u',
  intern = FALSE
)

# Fit the Heavy tailed model
objfun <- function(p) with(template_heavytail,list(value = fn(p),gradient = as.numeric(gr(p)),hessian = he(p)))
opt_heavy <- trust::trust(objfun,template_heavytail$par,2,50)
# coefficients
betas_heavy <- opt_heavy$argument[1:4]
alphas_heavy <- exp(opt_heavy$argument[5:6])

# random effects
template_heavytail_re$fn(opt_heavy$argument)
re_heavy <- matrix(with(template_heavytail_re$env,last.par.best[random]),ncol=2)


# predicted values
# point predictions- medians
state_point_predictions_heavy <- exp(betas_heavy[1] + re_heavy[ ,1])
plotdat_medians_heavy <- tibble(
  state = names(state_observations),
  obs = state_observations,
  pred_gaussian = state_point_predictions,
  pred_heavy = state_point_predictions_heavy
) %>%
  arrange(obs) %>%
  mutate(state = factor(state,levels = state))

median_plot_heavy <- plotdat_medians_heavy %>%
  ggplot(aes(x = state)) +
  theme_bw() +
  geom_hline(yintercept = median(plotdat_medians_heavy$obs)) +
  geom_point(aes(y = obs),size = 3) +
  geom_point(aes(y = pred_gaussian),pch = 4,size = 5) +
  geom_point(aes(y = pred_heavy),pch = 3,size = 5) +
  labs(x = "State sorted by median premium",y = "Median monthly premium ($)") +
  theme(text = element_text(size = PLOTTEXTSIZESMALL),axis.text.x = element_text(angle = 90))

ggsave(file.path(plotpath,"medians.pdf"),plot = median_plot_heavy,width = 7,height = 3.5)
ggsave(file.path(plotpath,"medians.tiff"),plot = median_plot_heavy,width = 7,height = 3.5,dpi=300)

# errors
mse_gaussian <- with(plotdat_medians_heavy,mean((obs-pred_gaussian)^2))
mse_t <- with(plotdat_medians_heavy,mean((obs-pred_heavy)^2))

(1 - mse_t/mse_gaussian)



