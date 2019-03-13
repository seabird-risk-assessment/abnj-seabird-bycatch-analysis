## source('../../../../../functions.r')
## includemk('../vars.mk')
## load(MAINDATA)
library(data.table)

source('../functions.r')
source('pst_parameter_chosen_values.r')
source('../config.r')

dat <- fread('generated/@estimates_table_final.csv')
dat <- dat[sp != 'DAQ']

spp <- fread('../input-data/acap-species-list.csv')
spp <- spp[base_species == T]
spp <- spp[code.sra != 'DAQ']

load('../input-data/species-demography/pop-corrections.rdata', verb=T)
setnames(smpls, 'sp', 'spcode')

d <- setdiff(spp$code.sra, dat$sp)
stopifnot(length(setdiff(spp$code.sra, dat$sp)) == 0)


## * Population size

pops <- readRDS('../input-data/species-demography/colony-population-sizes-samples.rds')
sp_pop_samples <- pops[, .(breeding_pairs = sum(breeding_pairs)), .(species_code, sample)]


stopifnot(length(c(setdiff(spp$code.sra, sp_pop_samples$species_code),
                   setdiff(sp_pop_samples$species_code, spp$code.sra))) == 0)

## * Samples

summt <- list()
ntot_summ <- list()

dir.create('generated/parameter_samples', showWarnings=F)

demsamples <- NULL

cat('\n=== Creating samples for each species\n')

i=2
for (i in seq_len(nrow(dat))) {
    d <- dat[i]
    sp <- d$sp
    cat(sp, '\t')
    params <- c('A','Scurr','Sopt','PB') 
    rm(list=params[params %in% ls()])

    ## A
    if (!is.na(d$A.mean) & !is.na(d$A.sd)) {
        A <- exp(rnorm(NSAMPLES, mean=log(d$A.mean), sd=d$A.sd)) 
        A_fixed <- d$A.mean
        summt[[sp]][['A']] <- sprintf('lognormal with mean %0.2f and sd %0.2f', log(d$A.mean), d$A.sd) 
    } else if (!is.na(d$A.min) & !is.na(d$A.max)) {
        A <- runif(NSAMPLES, min=d$A.min, max=d$A.max)
        A_fixed <- mean(c(d$A.min, d$A.max))
        summt[[sp]][['A']] <- sprintf('uniform btw %0.2f and %0.2f', d$A.min, d$A.max)
    }
    ##  Scurr
    if (!is.na(d$Scurr.mean) & !is.na(d$Scurr.sd)) {
        Scurr <- surv_var_norm_mean_se(NSAMPLES, smean=d$Scurr.mean, sse=d$Scurr.sd) 
        Scurr_fixed <- d$Scurr.mean
        summt[[sp]][['Scurr']] <- sprintf('normal on logit scale with mean %0.2f and sd %0.2f', d$Scurr.mean, d$Scurr.sd) 
    } else if (!is.na(d$Scurr.min) & !is.na(d$Scurr.max)) {
        Scurr <- runif(NSAMPLES, min=d$Scurr.min, max=d$Scurr.max)
        Scurr_fixed <- mean(c(d$Scurr.min, d$Scurr.max))
        summt[[sp]][['Scurr']] <- sprintf('uniform btw %0.2f and %0.2f', d$Scurr.min, d$Scurr.max)
    }

    ##  Sopt
    if (!is.na(d$Sopt.mean) & !is.na(d$Sopt.sd)) {
        Sopt <- surv_var_norm_mean_se(NSAMPLES, smean=d$Sopt.mean, sse=d$Sopt.sd) 
        Sopt_fixed <- d$Sopt.mean
        summt[[sp]][['Sopt']] <- sprintf('normal on logit scale with mean %0.2f and sd %0.2f', d$Sopt.mean, d$Sopt.sd) 
    } else if (!is.na(d$Sopt.min) & !is.na(d$Sopt.max)) {
        Sopt <- runif(NSAMPLES, min=d$Sopt.min, max=d$Sopt.max)
        Sopt_fixed <- mean(c(d$Sopt.min, d$Sopt.max))
        summt[[sp]][['Sopt']] <- sprintf('uniform btw %0.2f and %0.2f', d$Sopt.min, d$Sopt.max)
    }

    ##  PB
    if (!is.na(d$PB.mean) & !is.na(d$PB.sd)) {
        if (d$PB.mean != 1)
            PB <- surv_var_norm_mean_se(NSAMPLES, smean=d$PB.mean, sse=d$PB.sd) else
                                                                                       PB <- rep(d$PB.mean, NSAMPLES)
        PB_fixed <- d$PB.mean
        summt[[sp]][['PB']] <- sprintf('normal on logit scale with mean %0.2f and sd %0.2f', d$PB.mean, d$PB.sd) 
    } else if (!is.na(d$PB.min) & !is.na(d$PB.max)) {
        PB <- runif(NSAMPLES, min=d$PB.min, max=d$PB.max)
        PB_fixed <- mean(c(d$PB.min, d$PB.max))
        summt[[sp]][['PB']] <- sprintf('uniform btw %0.2f and %0.2f', d$PB.min, d$PB.max)
    }

    ##  N
    nads <- NA
    NN <- NN_nbp_fixed <- NN_a_fixed <- NN_s_fixed <- NN_pb_fixed <- list(all=NA, min=NA, s=NA, a=NA, r=NA, m=NA)
    breedpairs_tot <- sp_pop_samples[species_code == sp, breeding_pairs]
    summt[[sp]][['Nbp']] <- sprintf('from colonies')

    ## take lower quartile
    breedpairs_min <- breedpairs_tot[breedpairs_tot<=quantile(breedpairs_tot, 0.25)]
    if (length(breedpairs_min) < NSAMPLES) {
        breedpairs_min <- sample(breedpairs_min, NSAMPLES, replace=T)
    } else if (length(breedpairs_min) > NSAMPLES)
        breedpairs_min <- sample(breedpairs_min, NSAMPLES)
    breedpairs_tot <- sample(breedpairs_tot, NSAMPLES)
    totpairs_tot <- breedpairs_tot / PB
    breedpairs <- breedpairs_min ##
    breedpairs_fixed <- mean(breedpairs) 
    totpairs <- breedpairs / PB
    totpairs_nbp_fixed <- breedpairs_fixed / PB
    totpairs_pb_fixed <- breedpairs / PB_fixed

    totads <- totpairs * 2
    totads_tot <- totpairs_tot * 2
    totads_nbp_fixed <- totpairs_nbp_fixed * 2
    totads_pb_fixed <- totpairs_pb_fixed * 2
    NN     <- ntot2(Nadtot=totads,     surv=Scurr, afr=A, nsims=NSAMPLES)
    NN_tot <- ntot2(Nadtot=totads_tot, surv=Scurr, afr=A, nsims=NSAMPLES)
    NN_nbp_fixed <- ntot2(Nadtot=totads_nbp_fixed, surv=Scurr, afr=A, nsims=NSAMPLES)
    NN_pb_fixed <- ntot2(Nadtot=totads_pb_fixed, surv=Scurr, afr=A, nsims=NSAMPLES)
    NN_a_fixed <- ntot2(Nadtot=totads, surv=Scurr, afr=A_fixed, nsims=NSAMPLES)
    NN_s_fixed <- ntot2(Nadtot=totads, surv=Scurr_fixed, afr=A, nsims=NSAMPLES)

    ## ** Apply correction of Ntot
    corr <- smpls[spcode == spp[code.sra == sp, rho_group], ng_corr]
    NN_tot$all <- NN_tot$all * sample(corr, length(NN_tot$all), replace=T)

    df <- data.frame(sp=dat$sp[i], scurr=Scurr, sopt=Sopt, scurr_fixed=Scurr_fixed, sopt_fixed=Sopt_fixed,
                    a=A, a_fixed=A_fixed, pb=PB, pb_fixed=PB_fixed,
                    nbp=breedpairs, nbp_tot=breedpairs_tot, nbp_fixed=breedpairs_fixed, nads=totads,
                    nads_tot = totads_tot,
                    nads_nbp_fixed=totads_nbp_fixed, nads_pb_fixed=totads_pb_fixed,
                    m=NN$m, m_s_fixed=NN_s_fixed$m,
                    r=NN$r, r_a_fixed=NN_a_fixed$r, r_s_fixed=NN_s_fixed$r,
                    ntot=NN$all, ntot_tot=NN_tot$all, ntot_nbp_fixed=NN_nbp_fixed$all,
                    ntot_pb_fixed=NN_pb_fixed$all,
                    ntot_a_fixed=NN_a_fixed$all, ntot_s_fixed=NN_s_fixed$all)

    demsamples[[sp]] <- df

    ntot_summ[[sp]][['nbptot_p25']] <- quantile(df$nbp_tot, .25, na.rm=T)
    ntot_summ[[sp]][['nbptot_med']] <- quantile(df$nbp_tot, .5, na.rm=T)
    ntot_summ[[sp]][['nbptot_mean']] <- mean(df$nbp_tot, na.rm=T)
    ntot_summ[[sp]][['nbptot_lcl']] <- quantile(df$nbp_tot, .025, na.rm=T)
    ntot_summ[[sp]][['nbptot_ucl']] <- quantile(df$nbp_tot, .975, na.rm=T)
    
}

cat('\n\n')



save(demsamples, file='generated/dem-samples.rdata')
spl <- names(summt)
summt <- cbind(sp = spl, data.table(do.call('rbind', summt)))
ntot_summ <- cbind(sp = spl, data.table(do.call('rbind', ntot_summ)))
print(summt)
fwrite(summt, 'generated/samples_creation_summary.csv')
fwrite(ntot_summ, 'generated/ntot_summary_before_truncation.csv')

