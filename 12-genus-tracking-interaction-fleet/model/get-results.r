library(coda)
library(data.table)
library(parallel)
library(stringr)

options(scipen = 100)

load('../generated/model/model-data-image.rdata', v=F)
load('../generated/model/model-fitted.rdata', v=F)


## * Turn samples to long form
mcmc <- rbindlist(lapply(1:length(codasamples), function(ch) {
    d <- data.table(codasamples[[ch]])
    d[, chain := ch]
    d[, ch_sample := 1:.N]
    return(melt(d, id.vars = c('chain', 'ch_sample')))
}))

setorder(mcmc, variable, chain, ch_sample)

mcmc[, sample := 1L:.N, variable]

## * Get meaning out of the parameters

agg_fit[,  id := 1L:.N]
agg_pred[, id := 1L:.N]

get_meaning <- function(vars, covariate = NULL) {
    ## ** Returns the covariates given the variable name
    ## ** If `covariate` is NULL, returns all the covariates as data.table, otherwise returns a vector of single given covariate

    unvars <- unique(vars)
    
    if (!is.null(covariate) & length(covariate) > 1)  stop('get_meaning only accepts a single covariate to be specified')
    maxdim <- max(str_count(unvars, ',')) + 1L
    vartype <- sub('\\[.*\\]', '', unvars)
    ind1 <- rep(NA_integer_, length.out = length(unvars))
    singleidx <- grepl('.*\\[([0-9]+)\\].*', unvars)
    ind1[singleidx] <- as.integer(sub('.*\\[([0-9]+)\\].*', '\\1', unvars[singleidx]))
    ind2 <- rep(NA_integer_, length.out = length(unvars))
    if (maxdim > 1) {
        firstidx  <- grepl('\\[[0-9]+,', unvars)
        ind1[firstidx] <- as.integer(sub('.*\\[([0-9]+),.*', '\\1', unvars[firstidx]))
        secondidx <- grepl('\\[[0-9]+,[0-9]+\\]', unvars)
        ind2[secondidx] <- as.integer(sub('.*,([0-9]+)\\].*', '\\1', unvars[secondidx]))
    }
    dt <- data.table(unvars, vartype, ind1, ind2)
    
    add_attributes <- function(vartypes, from_dat, cols=NULL) {
        if (any(vartypes %in% dt$vartype)) {
            max1 <- sapply(vartypes, function(vt) {
                if (vt %in% dt$vartype) {
                    return(dt[vartype %in% vt, max(ind1)])
                } else return(0)
            })
            if (length(unique(setdiff(max1, 0))) != 1) {
                print(max1)
                stop('Inconsistencies between vartypes')
            }
            max1 <- unique(setdiff(max1, 0))
            if ('data.table' %in% class(from_dat)) {
                if (is.null(cols)) {
                    cols <- setdiff(names(from_dat), 'id')
                }
                stopifnot(max1 == max(from_dat$id))
                for (v in cols) {
                    dt[vartype %in% vartypes, eval(v) := from_dat[match(ind1, id), get(v)]]
                }
            } else {
                stopifnot(max1 == length(from_dat))
                dt[vartype %in% vartypes, eval(cols) := from_dat[ind1]]
            }
        }
    }
    print('Captures...')
    add_attributes(c('captures_o', 'dead_captures_o', 'live_captures_o', 'loglik_dead_o', 'loglik_live_o', 'q_g_o'),
                   agg_fit,
                   c('fishery_group', 'flag', 'species_group', 'species_code', 'site', 'grid_id', 'quarter'))
    
    print('APF')
    add_attributes(c('apf_t', 'observable_captures_t', 'observable_captures_ident_t', 'incidents_t'),
                   agg_pred,
                   c('fishery_group', 'flag', 'species_group', 'species_code', 'role', 'site', 'grid_id', 'quarter'))

    print('identifications')
    add_attributes(c('dead_unident_o', 'live_unident_o', 'loglik_unident_dead', 'loglik_unident_live',
                     'p_live_cap', 'q_f', 'q_f0'),
                   obs_fg, 'fishery_group')

    print('N')
    add_attributes(c('ntot', 'ntot_prior'), allspp, 'species_code')

    print('Vulnerability')
    add_attributes(c('q_g', 'q_g0', 'p_survive_cap'), attr(DUMMYVAR_SPECIES_GROUP_O, 'levels'))

    #print('Interaction')
    #add_attributes(c('q_gf', 'q_gf0'), attr(DUMMYVAR_INTERACTION_O, 'levels'))

    if (dt[grep('\\[', unvars), any(apply(.SD, 1, function(x) all(is.na(x)))), .SDcols = c('ind1', 'ind2')])
        stop('Some indices not parsed properly')

    nocovs <- dt[, apply(.SD, 1, function(x) all(is.na(x))), .SDcols = setdiff(names(dt), c('unvars', 'vartype', 'ind1', 'ind2'))]
    
    if (any(nocovs & dt[, !(vartype %in% c('q0', 'p_identified', 'p_observable', 'lp__', grep('^log_lik', vartype, val=T)))])) {
        warning('Some variables with no associated covariates')
        print( dt[nocovs==T & !(vartype %in% c('q0', 'p_identified', 'p_observable', 'lp__', grep('^log_lik', vartype, val=T))), .N, vartype] )
    }
    
    dt <- dt[match(vars, unvars)]
    setnames(dt, 'unvars', 'variable')
    
    if (!is.null(covariate)) {
        return(dt[, get(covariate)])
    } else return(dt)
    
}


vars <- levels(mcmc$variable)
vars <- vars[!grepl('_s$', sub('\\[.*\\]', '', vars))]

mc_attributes <- get_meaning(vars)
mc_attributes[grid_id %in% 0, grid_id := NA]
mc_attributes[quarter   %in% 0, quarter   := NA]

mcmc_core <- mcmc[grep('^[apf|capt]', variable, invert=TRUE)]
mcsumm <- mcmc[, .(mean = mean(value),
                  min  = min(value),
                  max  = max(value),
                  lcl  = quantile(value, 0.025, names=F),
                  ucl  = quantile(value, 0.975, names=F)),
              .(variable)]
mcsumm <- merge(mcsumm, mc_attributes, on = 'variable')

save(mcmc, mcmc_core, mc_attributes, mcsumm, file='../generated/model/mcmc-results.rdata', compress=F)
saveRDS(mcsumm, file='../generated/model/mcmc-results-summary.rds', compress=F)

