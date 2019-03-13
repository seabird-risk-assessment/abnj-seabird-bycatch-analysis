library(compiler)
library(data.table)

source('../../functions.r')

pst_rhos <- readRDS('../../input-data/species-demography/nu.rds')

load('../../input-data/species-demography/pop-corrections.rdata', verb=T)

load('../generated/dem-samples.rdata', verb=T)

sppchar <- fread('../acap-species-list.csv')
setkeyv(sppchar, 'code.sra')


pstfun <- function(sp1) {
    d <- demsamples[[sp1]]

    sp0 <- sp1

    A <- d$a
    S <- d$scurr
    nads <- 2 * d$nbp_tot / d$pb

    ## N/Nbp
    N_g <- round(nads * S^(1-A))
    
    ## Rmax from Niel&Lebreton
    rmax_nl <- Rmax_NL(d$sopt, d$a)

    ## PST with f=1
    pst_nlg <- 0.5 * rmax_nl * N_g

    attr(pst_nlg, 'rmax') <- rmax_nl
    attr(pst_nlg, 'N') <- N_g
    return(pst_nlg)
}

## with f=1
pst1 <- lapply(names(demsamples), pstfun)

pst1nlg <- as.data.frame(pst1)
names(pst1nlg) <- names(demsamples)



## * Add overall rho correction

rho <- 0.5
pst1nlgr <- pst1nlg * rho
rhovals <- rho


## * Turn to long form

pst_samples <- melt(as.data.table(pst1nlgr), measure.vars=names(pst1nlgr), variable.name='species_code', value.name = 'pst')
pst_samples[, sample := 1:.N, species_code]


## * Save
save(pst1nlg, pst1nlgr, rhovals, pst_samples, file = '../generated/pst-samples.rdata')


pstsumm <- rbindlist(lapply(names(pst1nlgr), function(sp) {
    x <- pst1nlgr[[sp]]
    data.table(sp=sp, mean = mean(x),
               lcl = quantile(x, 0.025, names = F),
               ucl = quantile(x, 0.975, names = F))
}))
save(pstsumm, file = '../generated/pst-summary.rdata')
