library(data.table)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

nchains <- 2

source('../../config.r')

## * Data
model_data <- new.env()
load('../generated/model/data.rdata', envir = model_data)
## attach(model_data)

## Different BETA0 between JAGS and Stan
model_data$BETA0 <- 1 #model_data$BETA0^-0.

## ## Overwrite the interaction variables, so the interaction is by species
## model_data$DUMMYVAR_INTERACTION_O <- model.matrix(data=data.frame(y=model_data$OVERLAP_O, 
##     x1=as.factor(model_data$SPECIES_O), 
##     x2=as.factor(model_data$FISHERY_GROUP_O)), y~x1:x2-1)

## model_data$N_INTERACTION_CAT_O <- ncol(model_data$DUMMYVAR_INTERACTION_O)

## model_data$DUMMYVAR_INTERACTION_T <- model.matrix(data=data.frame(y=model_data$OVERLAP_T, 
##     x1=as.factor(model_data$SPECIES_T), 
##     x2=as.factor(model_data$FISHERY_GROUP_T)), y~x1:x2-1)

## model_data$N_INTERACTION_CAT_T <- ncol(model_data$DUMMYVAR_INTERACTION_T)

save(model_data, file='../generated/model/model_data.rdata', compress = F)

warmup <- 2000
nsamples <- 1500

## * Run model
fit <- stan(file       = 'model.stan', 
           model_name = "SRAmodel",
           pars       = c('apf_t', 'q0', 'q_f0', 'q_g0', 'q_gf0',
                          'p_identified_bird',
                          'p_identified_family',
                          'captures_o',
			  'sigma_f',
			  'sigma_g',
                          'q_prior'),
           include    = TRUE, ## TRUE to only store parameters in `pars`
           data       = model_data, 
           iter       = warmup + nsamples,
           warmup     = warmup,
           control    = list(max_treedepth = 15, adapt_delta=0.8),
           chains     = nchains,
           diagnostic_file='../generated/model/diag.txt',
           verbose    = F,
           seed       = 3859285)


codasamples <- As.mcmc.list(fit)
mcmc <- data.table(do.call('rbind', codasamples))

save(codasamples, fit,
     file = '../generated/model/model-fitted.rdata', compress = F)

