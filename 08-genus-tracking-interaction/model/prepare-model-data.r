library(data.table)
library(sp)
library(sf)
library(countrycode)


source('../../functions.r', chdir=T)
source('../../config.r', chdir=T)

outdir <- '../generated/model'
dir.create(outdir, showWarnings = F)
           
load('../../input-data/cryptic-mortality/p_observable.rdata', v=T)

all_taxa <- fread('../acap-species-list.csv')
setnames(all_taxa, 'code.sra', 'species_code')

load('../generated/overlaps-and-captures.rdata', v=T)

print('Overlaps ... ')
overlap_obs[!is.na(countrycode(flag, 'iso3c', 'country.name', warn=F)), flag_name := countrycode(flag, 'iso3c', 'country.name', warn=F)]
overlap_tot[!is.na(countrycode(flag, 'iso3c', 'country.name', warn=F)), flag_name := countrycode(flag, 'iso3c', 'country.name', warn=F)]
effort_tot[!is.na(countrycode(flag, 'iso3c', 'country.name', warn=F)), flag_name := countrycode(flag, 'iso3c', 'country.name', warn=F)]

overlap_obs[flag %in% c('UNK', 'XX'), flag_name := 'Unknown']
overlap_tot[flag %in% c('UNK', 'XX'), flag_name := 'Unknown']
effort_tot[flag %in% c('UNK', 'XX'), flag_name := 'Unknown']

stopifnot(!any(is.na(overlap_obs$flag_name)) & !any(is.na(overlap_tot$flag_name)) & !any(is.na(effort_tot$flag_name)))



## * Merge overlaps

overlap_obs[, type := 'observed']
overlap_tot[, type := 'total']
# restrict to core area for prediction
overlap_tot <- overlap_tot[latitude < -20 & latitude > -50]

setnames(overlap_obs, 'observed_hooks', 'hooks')

stopifnot(all(setdiff(names(overlap_obs), names(overlap_tot)) == 'captures') & all(setdiff(names(overlap_tot), names(overlap_obs)) == 'site'))


## * Restrict years

overlap_obs  <- overlap_obs[year %in% YEARS_FIT]
overlap_tot  <- overlap_tot[year %in% YEARS_PRED]
effort_tot  <- effort_tot[yy %in% YEARS_PRED]

overlaps <- rbind(overlap_obs, overlap_tot, fill=T)

print('Species groups ... ')
## * Species groups

overlaps[all_taxa,
         `:=`(species_group = i.vulgroup,
	      family        = i.family,
              base_species  = i.base_species,
              taxa          = i.Taxa),
         on = 'species_code']

#stopifnot(nrow(overlaps[is.na(species_group)]) == 0)


print('fishery groups ... ')
## * Fishery groups

overlaps[,   fishery_group := flag]
effort_tot[, fishery_group := flag]


## * Check that all fishery groups for predictions are among fishery groups for vulnerability estimation
#if (length(setdiff(overlaps[type == 'total', fishery_group], overlaps[type == 'observed', fishery_group]))) {
#    stop('Some fishery groups in ov_pred not in ov_fit: ',
#         paste(setdiff(overlaps[type == 'total', fishery_group], overlaps[type == 'observed', fishery_group]), collapse=', '))
#}


## ** Keep track of fishery groups
#fgroups <- unique(overlaps[, .(fishery_group, flag, fishery)])
fgroups <-  unique(overlaps[, .(flag)])
fgroups[, slug := slugify(flag)]
fwrite(fgroups, '../generated/fishery-groups.csv')

overlaps0 <- copy(overlaps)


print('Remove rows with no overlap ... ')
## * Remove rows with no effort or no overlap
#overlaps <- overlaps[hooks > 0 & overlap > 0]
#overlaps[type=='observed', overlap := pmax(overlap, 1)]

overlaps[, overlap := overlap / 1e6]

overlap_unident <- overlaps[base_species == F]
overlap_ident   <- overlaps[base_species == T]
stopifnot(nrow(overlap_ident[is.na(species_group)]) == 0)

print('Make factors ... ')
## * Turn species, species group and fishery group to factor

allspp <- sort(unique(overlap_ident$species_code))
overlap_ident[, species_code := factor(species_code, levels = allspp)]
if (any(is.na(overlap_ident$species_code)))
    stop('Pb with species factor levels')

allspg <- sort(unique(overlap_ident$species_group)) # Species groups
overlap_ident[, species_group := factor(species_group, levels = allspg)]
if (any(is.na(overlap_ident$species_group)))
    stop('Pb with species group factor levels')

total_fg <- unique(overlaps$flag)
ident_fg <- unique(overlap_ident[type=='observed' & captures > 0, flag])
stopifnot(length(ident_fg) == 7)
unident_fg <- setdiff(unique(overlap_unident[type=='observed', flag]), ident_fg)
obs_fg <- sort(c(ident_fg, unident_fg))
allfg <- c(obs_fg, sort(setdiff(total_fg, obs_fg)))

overlap_ident[, fishery_group := factor(flag, levels = allfg)]
if (any(is.na(overlap_ident$fishery_group)))
    stop('Pb with fishery group factor levels')
overlap_unident[, fishery_group := factor(flag, levels = allfg)]
if (any(is.na(overlap_unident$fishery_group)))
    stop('Pb unidentified with fishery group factor levels')

allfamily <- rev(sort(unique(overlap_ident[!is.na(family) & nchar(family) > 0, family])))
overlap_ident[, family := factor(family, levels = allfamily)]
overlap_unident[, family := factor(family, levels = allfamily)]

ov_aves  <- overlap_unident[type == 'observed' & species_code == 'AVES']
ov_family <- overlap_unident[type == 'observed' & !is.na(family)]
ov_fit  <- overlap_ident[type == 'observed']
ov_pred <- overlap_ident[type == 'total']


print('Aggregate data ... ')
## * Aggregate data

## ** For estimation
fit_grid <- as.data.table(expand.grid(species_code=allspp, fishery_group=c(ident_fg, unident_fg)))
fit_grid[, species_code := factor(species_code, levels = allspp)]
fit_grid[, fishery_group := factor(fishery_group, levels = allfg)]
fit_grid <- merge(fit_grid, all_taxa[, .(species_code, species_group=vulgroup, family)], by='species_code')
fit_grid[, species_group := factor(species_group, levels = allspg)]
fit_grid[, family := factor(family, levels = allfamily)]

# Refine when you change the fishery_group !!!
fit_grid[, flag := fishery_group]

agg_fit <- ov_fit[, .(captures = sum(captures),
                      effort   = sum(hooks),
                      overlap  = sum(overlap),
                      nyears   = length(unique(year)),
                      grid_id  = 0L,
                      site     = 'All sites',
                      quarter    = 0L),
                  .(species_code, fishery_group)]

agg_fit <- merge(fit_grid, agg_fit, by=c('species_code', 'fishery_group'), all.x=TRUE)
## * Remove rows with no overlap
print('remove rows with tiny or zero overlap... ')
agg_fit <- agg_fit[overlap > 1e-6]
setorder(agg_fit, fishery_group, family, species_group, species_code)

agg_family <- ov_family[, .(captures = sum(captures)), .(fishery_group, family)]
setorder(agg_family, fishery_group, family)

agg_aves <- ov_aves[, .(captures = sum(captures)),
                  .(fishery_group)]
setorder(agg_aves, fishery_group)


print('Aggregate data for predictions... ')
## ** For predictions
agg_pred <- rbind(
    ## *** to compare predictions vs. observations (sum over years)
    cbind(agg_fit, role = 'obs_vs_pred'),
    ## *** spatial, all sites, no quarter (annual mean)
    ov_pred[latitude < -20 & latitude > -50,
            .(role = 'quarter0_spatial1', site = 'All sites', quarter=0L, overlap = sum(overlap) / length(YEARS_PRED)),
            .(species_code, species_group, fishery_group, grid_id)],
    ## *** spatial, by site, no quarter (annual mean)
    ov_pred[latitude < -20 & latitude > -50,
            .(role = 'quarter0_spatial2', quarter=0L, overlap = sum(overlap) / length(YEARS_PRED)),
            .(species_code, species_group, fishery_group, grid_id, site)],
    ## *** flag, not spatial (annual mean)
    ov_pred[latitude < -20 & latitude > -50,
            .(role = 'flag_quarter0_spatial0', site = 'All sites', grid_id=0L, quarter=0L, overlap = sum(overlap) / length(YEARS_PRED)),
            .(species_code, species_group, fishery_group, flag)],
    fill=T)

## * Remove rows with no overlap
print('remove rows with no overlap... ')
agg_pred <- agg_pred[overlap > 0 & is.finite(overlap)]

## * Model inputs

## ** Species and groups

print('Save for model... ')
N_SPECIES  <- length(allspp)
N_SPECIES_VUL_GROUP <- length(allspg)
N_FISHERY_GROUP <- length(allfg)
N_OBS_FISHERY_GROUP <- max(as.integer(agg_fit$fishery_group))
stopifnot(N_OBS_FISHERY_GROUP == length(obs_fg))
N_FAMILY_GROUP <- length(allfamily)


## ** P_observable
#betapars <- estBetaParams(mean(p_obs), var(p_obs))
#POBS_BETA_ALPHA <- betapars$alpha
#POBS_BETA_BETA  <- betapars$beta


## *** Sort data and index groups to sum over in model
setorder(agg_fit, fishery_group, family, species_group, species_code)
## ** Unidentified

idx <- agg_fit[, .(start=min(.I), end=max(.I)), .(fishery_group)]
idx <- merge(idx, agg_aves, by = 'fishery_group', all.x=T, all.y=F)

setorder(idx, fishery_group)
AVES <- agg_aves$captures
AVES_N <- nrow(idx)
stopifnot(AVES_N == N_OBS_FISHERY_GROUP)
AVES_FISHERY_O <- as.numeric(idx$fishery_group)
AVES_START_O <- idx$start
AVES_END_O   <- idx$end
idx[, aves_prior := 10000]
idx[fishery_group %in% ident_fg, aves_prior := 1]
AVES_PRIOR <- idx$aves_prior

idx <- agg_fit[, .(start=min(.I), end=max(.I)), .(fishery_group, family)][order(fishery_group, family)]
idx <- merge(idx, agg_family, by=c('fishery_group', 'family'))
FAMILY <- idx$captures
FAMILY_FISHERY_N <- nrow(idx)
stopifnot(length(FAMILY)==FAMILY_FISHERY_N)
FAMILY_FISHERY_O <- as.numeric(idx$fishery_group)
FAMILY_FAMILY_O <- as.numeric(idx$family)
FAMILY_START_O <- idx$start
FAMILY_END_O   <- idx$end


## *** Observed for estimation
N_ROW_O <- nrow(agg_fit)
OVERLAP_O <- agg_fit[, overlap]
CAPTURES_O <- agg_fit[, captures]
SPECIES_O <- agg_fit[, as.integer(species_code)]
SPECIES_GROUP_O <- agg_fit[, as.integer(species_group)]
FISHERY_GROUP_O <- agg_fit[, as.integer(fishery_group)]
FAMILY_GROUP_O <- agg_fit[, as.integer(family)]
Q_SCALE_SPECIES_GROUP <- agg_fit[, .(scale=sum(captures)/sum(overlap)), species_group]$scale
Q_SCALE_O <- Q_SCALE_SPECIES_GROUP[SPECIES_GROUP_O]

## *** Total for prediction
N_ROW_T <- nrow(agg_pred)
OVERLAP_T <- agg_pred[, overlap]
SPECIES_T <- agg_pred[, as.integer(species_code)]
SPECIES_GROUP_T <- as.integer(agg_pred[, species_group])
FISHERY_GROUP_T <- agg_pred[, as.integer(fishery_group)]
Q_SCALE_T <- Q_SCALE_SPECIES_GROUP[SPECIES_GROUP_T]


DUMMYVAR_SPECIES_GROUP_O <- model.matrix(data=data.frame(y=OVERLAP_O, x=as.factor(SPECIES_GROUP_O)), y~x-1)
attr(DUMMYVAR_SPECIES_GROUP_O, 'levels') <- data.table(id = 1:N_SPECIES_VUL_GROUP, species_group = allspg)

DUMMYVAR_FISHERY_GROUP_O <- model.matrix(data=data.frame(y=OVERLAP_O, x=as.factor(FISHERY_GROUP_O)), y~x-1)
attr(DUMMYVAR_FISHERY_GROUP_O, 'levels') <- data.table(id = 1:N_FISHERY_GROUP, fishery_group = allfg)

DUMMYVAR_INTERACTION_O <- model.matrix(data=data.frame(y=OVERLAP_O, x1=as.factor(SPECIES_GROUP_O), x2=as.factor(FISHERY_GROUP_O)), y~x1:x2-1)
attr(DUMMYVAR_INTERACTION_O, 'levels') <- data.table(id = 1:(N_SPECIES_VUL_GROUP*N_OBS_FISHERY_GROUP),
                                                     expand.grid(species_group=allspg, fishery_group=obs_fg, stringsAsFactors=F))
    
N_INTERACTION_CAT_O <- ncol(DUMMYVAR_INTERACTION_O)


DUMMYVAR_SPECIES_GROUP_T <- model.matrix(data=data.frame(y=OVERLAP_T, x=as.factor(SPECIES_GROUP_T)), y~x-1)
attr(DUMMYVAR_SPECIES_GROUP_T, 'levels') <- data.table(id = 1:N_SPECIES_VUL_GROUP, species_group = allspg)

DUMMYVAR_FISHERY_GROUP_T <- model.matrix(data=data.frame(y=OVERLAP_T, x=as.factor(FISHERY_GROUP_T)), y~x-1)
attr(DUMMYVAR_FISHERY_GROUP_T, 'levels') <- data.table(id = 1:N_FISHERY_GROUP, fishery_group = allfg)

DUMMYVAR_INTERACTION_T <- model.matrix(data=data.frame(y=OVERLAP_T, x1=as.factor(SPECIES_GROUP_T), x2=as.factor(FISHERY_GROUP_T)), y~x1:x2-1)
attr(DUMMYVAR_INTERACTION_T, 'levels') <- data.table(id = 1:(N_SPECIES_VUL_GROUP*N_FISHERY_GROUP),
                                                     expand.grid(species_group=allspg, fishery_group=allfg, stringsAsFactors=F))

N_INTERACTION_CAT_T <- ncol(DUMMYVAR_INTERACTION_T)
## * Save

save.image(file=file.path(outdir, 'model-data-image.rdata'), compress = F)

save(N_SPECIES, N_SPECIES_VUL_GROUP, N_FISHERY_GROUP, N_OBS_FISHERY_GROUP,
     CAPTURES_O, 
     N_ROW_O, N_ROW_T, N_FAMILY_GROUP,
     FISHERY_GROUP_O, FISHERY_GROUP_T, SPECIES_O, SPECIES_T,
     SPECIES_GROUP_O, OVERLAP_O, 
     SPECIES_GROUP_T, OVERLAP_T,
     Q_SCALE_O, Q_SCALE_T,
     AVES, AVES_N, AVES_FISHERY_O, AVES_START_O, AVES_END_O, AVES_PRIOR,
     FAMILY, FAMILY_FISHERY_N, FAMILY_GROUP_O, FAMILY_FISHERY_O, FAMILY_FAMILY_O, FAMILY_START_O, FAMILY_END_O,
     DUMMYVAR_SPECIES_GROUP_O, DUMMYVAR_FISHERY_GROUP_O, DUMMYVAR_INTERACTION_O,
     N_INTERACTION_CAT_O,
     N_INTERACTION_CAT_T,
     DUMMYVAR_SPECIES_GROUP_T, DUMMYVAR_FISHERY_GROUP_T, DUMMYVAR_INTERACTION_T,
     file = file.path(outdir, 'data.rdata'), compress = F)


