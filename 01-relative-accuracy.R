# Relative accuracy of the BCT approximations #


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
  'microbenchmark',
  'parallel'
)
for (pkg in pkgs) {
  if (!require(pkg,character.only = TRUE,quietly = TRUE)) {
    cat(paste0("Could not find package ",pkg,", installing from CRAN.\n"))
    install.packages(pkg)
    require(pkg,character.only = TRUE,quietly = TRUE)
  }
}

options(mc.cores = parallel::detectCores())

PLOTTEXTSIZE <- 18
plotpath <- "~/work/projects/bct-ad/figures"

# Three functions
# Without loss of generality, fix sigma = 1
bc <- function(nu,yovermu) {
  out <- numeric(length(nu))
  for (i in 1:length(nu)) {
    if (nu[i] == 0) {
      out[i] <- log(yovermu)
    } else {
      out[i] <- ( (yovermu)^nu[i] - 1 ) / (nu[i])
    }
  }
  out
}
bc_series <- function(nu,M,yovermu) {
  # M:number of terms
  out <- numeric(length(nu))
  for (i in 1:length(nu)) {
    out[i] <- log(yovermu) + sum(nu[i]^(1:(M-1)) * (log(yovermu))^(2:M) / factorial(2:M))
  }
  out
}
# GL quadrature
# Compute the nodes and weights only once
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

bc_gl3 <- function(nu,yovermu) (yovermu-1) * sum(( ( (yovermu-1)*nn3 + 1 )^(nu-1) )*ww3)
bc_gl9 <- function(nu,yovermu) (yovermu-1) * sum(( ( (yovermu-1)*nn9 + 1 )^(nu-1) )*ww9)
bc_gl15 <- function(nu,yovermu) (yovermu-1) * sum(( ( (yovermu-1)*nn15 + 1 )^(nu-1) )*ww15)
bc_gl25 <- function(nu,yovermu) (yovermu-1) * sum(( ( (yovermu-1)*nn25 + 1 )^(nu-1) )*ww25)

relerror <- function(true,approx) log(abs(1 - true/approx),base=10)
relerrorstable <- function(true,approx) log(abs(expm1(log(approx) - log(true))),base=10)

## Plots of relative error ##
paramstotry <- expand.grid(
  nu = seq(-3,3,by=.01),
  yovermu = seq(.1,2,by=.01)
)
bcout <- seriesout4 <- seriesout10 <- seriesout25 <- seriesout100 <- glout3 <- glout9 <- glout15 <- glout25 <- numeric(nrow(paramstotry))
for (i in 1:nrow(paramstotry)) {
  bcout[i] <- with(paramstotry[i, ,drop=FALSE],bc(nu,yovermu))
  seriesout4[i] <- with(paramstotry[i, ,drop=FALSE],bc_series(nu,4,yovermu))
  seriesout10[i] <- with(paramstotry[i, ,drop=FALSE],bc_series(nu,10,yovermu))
  seriesout25[i] <- with(paramstotry[i, ,drop=FALSE],bc_series(nu,25,yovermu))
  seriesout100[i] <- with(paramstotry[i, ,drop=FALSE],bc_series(nu,100,yovermu))
  glout3[i] <- with(paramstotry[i, ,drop=FALSE],bc_gl3(nu,yovermu))
  glout9[i] <- with(paramstotry[i, ,drop=FALSE],bc_gl9(nu,yovermu))
  glout15[i] <- with(paramstotry[i, ,drop=FALSE],bc_gl15(nu,yovermu))
  glout25[i] <- with(paramstotry[i, ,drop=FALSE],bc_gl25(nu,yovermu))
}

errormat <- paramstotry
errormat$bc <- bcout
errormat$series4 <- seriesout4
errormat$series10 <- seriesout10
errormat$series25 <- seriesout25
errormat$series100 <- seriesout100
errormat$gl3 <- glout3
errormat$gl9 <- glout9
errormat$gl15 <- glout15
errormat$gl25 <- glout25

errormat <- within(errormat,{
  rel_error_series4 = relerror(bc,series4)
  rel_error_series10 = relerror(bc,series10)
  rel_error_series25 = relerror(bc,series25)
  rel_error_series100 = relerror(bc,series100)
  
  rel_error_gl3 = relerror(bc,gl3)
  rel_error_gl9 = relerror(bc,gl9)
  rel_error_gl15 = relerror(bc,gl15)
  rel_error_gl25 = relerror(bc,gl25)
  
})
# maxima
errormat <- as_tibble(errormat)
maxmat <- errormat %>%
  group_by(yovermu) %>%
  summarize(across(-contains("nu"),max)) %>%
  pivot_longer(contains("rel_error"),names_to = "type",values_to = "log_rel_error") %>%
  mutate(
    type = str_extract(type,"[a-z]+[0-9]+"),
    method = str_extract(type,"[a-z]+"),
    tuning = str_extract(type,"[0-9]+")
  ) %>%
  dplyr::select(yovermu,type,method,tuning,log_rel_error)


glplot <- maxmat %>%
  filter(method == 'gl') %>%
  ggplot(aes(x = yovermu,y = log_rel_error)) + 
  theme_bw() + 
  geom_line(aes(linetype = tuning)) +
  scale_linetype_manual(values = c("3" = "longdash","9" = "dashed","15" = "dotdash","25" = "solid")) +
  geom_hline(yintercept = log(.Machine$double.eps,base=10),linetype = 'dotted') +
  coord_cartesian(ylim = c(-16,5)) +
  theme(text = element_text(size = PLOTTEXTSIZE)) + 
  labs(x = expression(Y/mu),y = expression("log"["10"]*"(rel. error.)"),linetype = "k")


seriesplot <- maxmat %>%
  filter(method == 'series') %>%
  ggplot(aes(x = yovermu,y = log_rel_error)) + 
  theme_bw() + 
  geom_line(aes(linetype = tuning)) +
  scale_linetype_manual(values = c("4" = "longdash","10" = "dashed","25" = "dotdash","100" = "solid")) +
  geom_hline(yintercept = log(.Machine$double.eps,base=10),linetype = 'dotted')+
  coord_cartesian(ylim = c(-16,5)) +
  theme(text = element_text(size = PLOTTEXTSIZE)) + 
  labs(x = expression(Y/mu),y = expression("log"["10"]*"(rel. error.)"),linetype = "m")

ggsave(glplot,file = file.path(plotpath,"glplot.pdf"),width=7,height=7)
ggsave(seriesplot,file = file.path(plotpath,"seriesplot.pdf"),width=7,height=7)

## Computational times ##
get_relative_times <- function(nu,yovermu,times = 100) {
  times <- microbenchmark::microbenchmark(
    bc(nu,yovermu),
    bc_series(nu,4,yovermu),
    bc_series(nu,10,yovermu),
    bc_series(nu,25,yovermu),
    bc_series(nu,100,yovermu),
    bc_gl3(nu,yovermu),
    bc_gl9(nu,yovermu),
    bc_gl15(nu,yovermu),
    bc_gl25(nu,yovermu),
    times = times
  )
  times <- as_tibble(as.data.frame(times)) %>%
    mutate(
      method = str_extract(expr,"bc_*[a-z]*"),
      tuning = str_extract(expr,"[0-9]+")
    )
  times_bc <- times %>% filter(method == 'bc')
  times_approx <- times %>% filter(method != 'bc')
  bc_median <- median(times_bc$time)
  timesummary <- times_approx %>%
    group_by(method,tuning) %>%
    summarize(mediantime = median(time),.groups = "drop") %>%
    mutate(relative_median = 100 * mediantime / bc_median)
  
  bc_summary <- times_bc %>%
    group_by(method,tuning) %>%
    summarize(mediantime = median(time),.groups = "drop") %>%
    mutate(relative_median = 100 * mediantime / bc_median)
  
  timesummary %>%
    bind_rows(bc_summary)
  
}

paramstotry_time <- expand.grid(
  nu = seq(-3,3,by=.1),
  yovermu = seq(.1,2,by=.1)
)

timelist <- vector(mode = 'list',length = nrow(paramstotry_time))
for (i in 1:nrow(paramstotry_time)) timelist[[i]] <- paramstotry_time[i, ,drop=FALSE]
dotime <- function(lst) {
  out <- with(lst,get_relative_times(nu,yovermu))
  out$nu <- lst$nu
  out$yovermu <- lst$yovermu
  out
}
timesthatweretried <- mclapply(timelist,dotime)
timeframe <- bind_rows(timesthatweretried)
timeframe %>%
  group_by(method,tuning) %>%
  summarize(
    mean_median = mean(mediantime),
    sd_median = sd(mediantime),
    min_median = min(mediantime),
    max_median = max(mediantime),
    median_median = median(mediantime)
  )

# Profile over each
timeframe_profile_nu <- timeframe %>%
  group_by(method,tuning,yovermu) %>%
  summarize(mediantime = median(mediantime))

timeframe_profile_nu %>%
  mutate(type = paste(method,tuning,sep="_")) %>%
  ggplot(aes(x=yovermu,y=mediantime)) +
  theme_bw() +
  geom_line(aes(linetype = type))

timeframe_profile_yovermu <- timeframe %>%
  group_by(method,tuning,nu) %>%
  summarize(mediantime = median(mediantime))

timeframe_profile_yovermu %>%
  mutate(type = paste(method,tuning,sep="_")) %>%
  ggplot(aes(x=nu,y=mediantime)) +
  theme_bw() +
  geom_line(aes(linetype = type))

# No clear pattern or difference.
# Report results for yovermu = 1.5, nu = 1, based on 1MM runs
final_times <- get_relative_times(1,1.5,1e06)
# Relative to each other
final_times_gl <- final_times[1:4, ] %>% mutate(tuning = as.numeric(tuning)) %>% arrange(tuning)
final_times_series <- final_times[5:8, ] %>% mutate(tuning = as.numeric(tuning)) %>% arrange(tuning)

100 * final_times_series$mediantime / final_times_gl$mediantime





