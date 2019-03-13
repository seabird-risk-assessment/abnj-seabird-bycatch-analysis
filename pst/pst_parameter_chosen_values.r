## for survival
ssdgood <- 0.01
ssdmed <- 0.02
ssdpoor <- 0.03

## for age at first reproduction
a_mult_minonly <- expression(5/3)
a_mult_maxonly <- expression(3/5)
a_multmin_meanonly <- expression(3/4)
a_multmax_meanonly <- expression(5/4)

## for population size
nsdgood <- 0.1
nsdmed <- 0.2
nsdpoor <- 0.3
n_mult_minonly_max <- 3
n_mult_minonly_min <- 0.7
n_mult_maxonly_min <- 0.2
n_mult_maxonly_max <- 1.2

## for proportion of breeders
pb_sd_meanonly <- 0.05
PBan <- list(mean=0.9, sd=pb_sd_meanonly)
PBbi <- list(mean=0.6, sd=pb_sd_meanonly)
PBmix <- list(mean=0.75, sd=pb_sd_meanonly)


## number of bootstraps for sensitivity
nboots <- 500

