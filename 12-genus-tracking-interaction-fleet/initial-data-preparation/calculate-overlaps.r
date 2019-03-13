library(data.table)
library(sf)
library(countrycode)
source('../../config.r')

spp <- fread('../acap-species-list.csv')

load('../../input-data/grid/grids.rdata', v=T)

griddt <- as.data.table(grid_noland)
setnames(griddt, tolower(names(griddt)))

## * BIRD DENSITIES
alldensities <- fread('../generated/gridded-distribution-by-species.csv')


## * OBSERVED DATA
obs <- fread('../../input-data/effort-captures/combined_species.csv')
obs <- obs[source %in%  obsflag]
obs[, key := paste(source, row)]

# Remove taxa that are not in the five principal genera
outliers <- c('PUG', 'PFC', 'PQW', 'DAC', 'DAQ', 'DBN+DIW+DQS', 'DIX+DBN+DIW',
	'DPK', 'PUC', 'SVB', 'WANDERING') 


write.csv(obs[taxon %in% outliers], file='removed-records.csv')
obs <- obs[!(taxon %in% outliers)]

#recode
stopifnot(!(any(c('THA', 'DIO', 'PHX') %in% unique(obs$taxon))))
recode <- fread('recode.csv')
obs <- merge(obs, recode, by.x='taxon', by.y='code', all.x=T, all.y=F)
obs[!is.na(update), taxon := update]
stopifnot(all(obs$taxon %in% spp$code.sra)) 
obscaptures <- obs[, .(captures=sum(captures)), .(key, species_code=taxon)]

obseffort <- fread('../../input-data/effort-captures/combined_total.csv')
obseffort <- obseffort[source %in%  obsflag]
stopifnot(length(obsflag) == length(unique(obseffort$source)))
obseffort[, key := paste(source, row)]

# Now expand out to include  records for all species on all fishing effort
overlap_grid <- as.data.table(expand.grid(species_code=spp$code.sra, key=obseffort$key))
obsdata <- merge(overlap_grid, obseffort, by='key')
obsdata <- merge(obsdata, obscaptures, by=c('key', 'species_code'), all.x=T)
obsdata[is.na(captures), captures := 0]
print(obsdata[,(.N), species_code])

## ** Round to avoid floating point issues
obsdata[, `:=`(longitude = round(longitude, 1),
               latitude  = round(latitude, 1))]
alldensities[, `:=`(centre.x = round(centre.x, 1),
                    centre.y = round(centre.y, 1))]

stopifnot(obsdata[, all(longitude >= -180 & longitude <= 180 & latitude <= 0 & latitude >= -90)])

## ** Add grid_id
obsdata[as.data.table(grid_noland), grid_id := i.grid_id, on = c('longitude' = 'centre.x', 'latitude' = 'centre.y')]

if (obsdata[, any(is.na(grid_id))]) {
    miss <- obsdata[is.na(grid_id)]
    warning('Some observed data with no grid_id (on land?). Will be removed.')
    print(miss[, .(.N, eff = sum(observed_hooks)), .(longitude, latitude)])
    fwrite(miss, '../generated/observed-data-without-gridid.csv')
} else if (file.exists('../generated/observed-data-without-gridid.csv')) {
    file.remove('../generated/observed-data-without-gridid.csv')
}

## ** Get overlap

#overlap_obs:
# flag longitude latitude grid_id year quarter fishery observed_hooks captures density  overlap species_code


dens <- alldensities[, .(density = sum(density)), .(species_code=code.sra, quarter, centre.x, centre.y, grid_id, taxa)]
obsdata[dens, `:=`(density = i.density), on = c('species_code', 'grid_id', 'quarter')]
obsdata[is.na(density), density := 0]
obsdata[, overlap := observed_hooks * density]
overlap_obs <- obsdata[, .(flag=source, longitude, latitude, grid_id, year, quarter, fishery, observed_hooks, captures, density,  overlap, species_code)]

## ** Check captures with zero densities

captures_with_zero_density <- overlap_obs[captures > 0 & density == 0]

if (nrow(captures_with_zero_density)) {
    warning(captures_with_zero_density[, sum(captures)], ' captures with zero bird density (',
            round(100*captures_with_zero_density[, sum(captures)] / overlap_obs[, sum(captures)], 3), '% of all captures). Saved to `generated/captures-with-zero-density.csv`.')
    print(captures_with_zero_density[, sum(captures), species_code])
    fwrite(captures_with_zero_density, '../generated/captures-with-zero-density.csv')
}

effort_obs <- obsdata


## * TOTAL EFFORT

effort_tot_rfmo <- fread("../../input-data/effort-captures/effort-total-RFMO.csv",
                         na.strings = 'N/A') # because of Namibia iso-2 = NA

## ** Replace NZ effort
nz_effort <- fread('../../input-data/effort-captures/NZL/total-effort.csv')
effort_tot_rfmo <- rbind(effort_tot_rfmo[is.na(flag) | !(flag %in% 'NZ')], nz_effort)


## ** Round to avoid floating point issues
effort_tot_rfmo[, `:=`(longitude = round(lon5, 1),
                       latitude  = round(lat5, 1))]

effort_tot_rfmo[longitude > 180,  longitude := longitude - 360]
effort_tot_rfmo[longitude < -180, longitude := longitude + 360]

stopifnot(effort_tot_rfmo[, all(longitude >= -180 & longitude <= 180 & latitude <= 0 & latitude >= -90)])

#Standardise to three letter codes
effort_tot_rfmo[!is.na(countrycode(flag, 'iso2c', 'iso3c', warn=F)), flag := countrycode(flag, 'iso2c', 'iso3c', warn=F)]

## ** Add grid_id
effort_tot_rfmo[as.data.table(grid_noland), grid_id := i.grid_id, on = c('longitude' = 'centre.x', 'latitude' = 'centre.y')]

if (effort_tot_rfmo[, any(is.na(grid_id))]) {
    miss <- effort_tot_rfmo[is.na(grid_id)]
    warning('Some RFMO data with no grid_id (on land?). Will be removed.')
    effort_tot_rfmo <- effort_tot_rfmo[!is.na(grid_id)]
    print(miss[, .(.N, eff = sum(hooks)), .(longitude, latitude)])
    fwrite(miss, '../generated/rfmo-data-without-gridid.csv')
} else if (file.exists('../generated/rfmo-data-without-gridid.csv')) {
    file.remove('../generated/rfmo-data-without-gridid.csv')
}

# Check the merged overlap and density!!

## ** Get overlap and turn to long form
#overlap_tot:
# flag fishery longitude latitude grid_id year quarter  hooks  density  site    overlap species_code

alldensities[is.na(density), density := 0]
overlap_tot <- merge(alldensities[, .(grid_id, quarter, code.sra, taxa, site, centre.x, centre.y, density)], 
	effort_tot_rfmo[hooks > 0 & yy > 2011, .(grid_id, year=yy, quarter, flag, fishery, hooks)], 
	by=c('grid_id', 'quarter'), 
	allow.cartesian=TRUE)
overlap_tot[is.na(density), density := 0]
overlap_tot[, overlap := hooks * density]
overlap_tot <- overlap_tot[, .(flag, fishery, longitude=centre.x, latitude=centre.y, grid_id, year, quarter,  hooks,  density,  site,    overlap, species_code=code.sra)]

effort_tot <- effort_tot_rfmo[, -c('yq', 'lat5', 'lon5'), with=F]


save(overlap_obs, effort_obs, overlap_tot, effort_tot, file = '../generated/overlaps-and-captures.rdata', compress=F)
