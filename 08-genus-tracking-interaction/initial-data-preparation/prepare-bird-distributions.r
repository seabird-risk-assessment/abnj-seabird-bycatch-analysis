suppressPackageStartupMessages({
    library(data.table)
    library(raster)
    library(sf)
})

dir.create('generated', showWarnings=F, recursive = T)

sppchar <- fread('../acap-species-list.csv')[order(order)]
source('../../config.r')
source('../../functions.r')

load('../../input-data/grid/grids.rdata')

rhos <- readRDS('../generated/rhos.rds')

sp_cols <- read_sf('../../input-data/species-distributions/tracks_sp_col.gpkg')

## ** Get grid id
sp_cols <- st_wrap_dateline(sp_cols)
sp_cols <- st_join(sp_cols, grid)

ignore_species <- c("Waved Albatross", "Pink-footed Shearwater")
sp_cols <- as.data.table(sp_cols)
sp_cols <- sp_cols[!(species %in% ignore_species)]
sp_cols[species == 'Campbell Albatross', species := 'Campbell Black-browed Albatross']
sp_cols[species == 'Gibson\'s Albatross', species := 'Antipodean Albatross']
sp_cols[species == 'Light-mantled Albatross', species := 'Light-mantled Sooty Albatross']
sp_cols[species == 'Chatham Albatross', species := 'Chatham Island Albatross']

## * Population sizes

load('../generated/dem-samples.rdata', v=F)
dem <- rbindlist(lapply(names(demsamples), function(sp) {
    dem <- demsamples[[sp]]
    setDT(dem)
    dem[, sample := 1:.N]
    return(dem)
}))

pops <- readRDS('../../input-data/species-demography/colony-population-sizes-samples.rds')

pops[dem, `:=`(pb    = i.pb,
               scurr = i.scurr,
               afr     = i.a),
     on = c('species_code' = 'sp', 'sample')]


## ** Total population size using Gilbert's
pops[, ntot_uncorr := (2 * breeding_pairs / pb) * scurr ^ ( 1 - afr )]

## ** Apply correction of Ntot
pops[sppchar, rho_group := i.rho_group, on = c('species_code' = 'code.sra')]
pops[rhos, rho := i.rho, on = c('rho_group' = 'sptype', 'sample')]
pops[, ntot := ntot_uncorr * rho]

pops_mean <- pops[, .(ntot = round(mean(ntot))), .(species_code, site)]



## * Tracking distributions

dist <- readRDS('../../input-data/species-distributions/tracking-distributions-densities.rds')

sppchar[, cname2 := tolower(cname)]

lkup <- dist[, .(map_cname = unique(species))]
lkup[, map_cname2 := tolower(map_cname)]
lkup[map_cname2 == 'chatham albatross', map_cname2 := 'chatham island albatross']
lkup[map_cname2 == 'light-mantled albatross', map_cname2 := 'light-mantled sooty albatross']
lkup[sppchar, code.sra := i.code.sra, on = c('map_cname2' = 'cname2')]
stopifnot(!any(is.na(lkup$code.sra)))

dist[lkup, code.sra := i.code.sra, on = c('species' = 'map_cname')]
stopifnot(!any(is.na(dist$code.sra)))

dist[, source := 'Tracking distribution']


## * Use other species' distribution in some cases


## ** Copy single combination of species/site

site_prox <- function(dist, sp_ori, site_ori, sp_targ=sp_ori, site_targ) {
    d <- dist[code.sra == sp_ori & site == site_ori]
    d[, `:=`(code.sra = sp_targ,
             site = site_targ,
             source = 'Proxy colony distribution')]
    rbind(dist[!(code.sra == sp_targ & site == site_targ)], d)
}

## *** Atl. Yellow-nosed: Gough -> Tristan
dist <- site_prox(dist, sp_ori='DCR', site_ori='Gough Island', site_targ='Tristan da Cunha')


## * Species with no distribution at all -> use range maps

missing <- setdiff(pops_mean$species_code, dist$code.sra)

## ** Manual changes
missing <- c(
    ## *** Do not use range maps for these species
    setdiff(missing,
            'DIP'), # southern royal
    ## *** Force the use of range maps for these species
    c('MAI', 'MAH', # giant petrels
      'PHE', # light-mantled sooty albatross
      'PHU', # Sooty albatross
      'PRO', # white-chinned petrel
      'PCI', # grey petrel
      'DIB', # Buller's albatross
      'DIC'  # grey-headed albatross
      )
) 

if (length(missing)) {
    warning('Tracking distribution missing for:\n', sppchar[code.sra %in% missing, paste(cname, collapse=', ')],
            '\nRange maps will be used for them\n')

    rangemaps <- readRDS('../../input-data/species-distributions/birdlife-range-maps-densities.rds')
    rangemaps <- rangemaps[code.sra %in% missing]

    rangemaps <- rbindlist(lapply(1:4, function(q) {
        rangemaps[, .(code.sra, site = 'Range map', quarter = q, grid_id, centre.x, centre.y, density_norm,
                      area_water = range_area,
                      source = 'Range map')]
    }))

    dist <- rbind(dist[!(code.sra %in% missing),
                       .(code.sra, site, quarter, grid_id, centre.x = x, centre.y = y, density_norm, area_water, source)],
                  rangemaps, fill=T)
    pops_mean_missing <- pops_mean[species_code %in% missing]
    pops_mean_missing <- pops_mean_missing[, .(site = 'Range map', ntot = sum(ntot)), .(species_code)]
    pops_mean <- rbind(pops_mean[!(species_code %in% missing)],
                       pops_mean_missing, fill=T)
}


## * Create distribution for colonies with no distribution (but species has some colonies with distribution)

sp_cols[, cname2 := tolower(species)]
sp_cols[sppchar, code.sra := i.code.sra, on = 'cname2']
stopifnot(sp_cols[, sum(is.na(code.sra)) == 0])

dist_sp_site <- dist[, .(.N), .(species_code=code.sra, site)]

pops_mean[, p_world := ntot / sum(ntot), species_code]


## ** Missing species/colonies with no distribution

missing <- pops_mean[!dist_sp_site, on = c('species_code', 'site')]

## ** Distance between grid cell centers
griddt <- as.data.table(grid)
grid_pts <- st_as_sf(griddt[, -'geom', with=F], coords = c('centre.x', 'centre.y'), crs = 4326)
grid_pts_dist <- st_distance(grid_pts, grid_pts, by_element=F)

grid_nl <- as.data.table(grid_noland)


## ** Normalised density from each species/missing colony
dist_missing <- rbindlist(lapply(seq_len(nrow(missing)), function(i) {
    
    col <- sp_cols[missing[i], on = c('code.sra'='species_code', 'loc_new'='site')]
    dist_to_col <- as.numeric(grid_pts_dist[col$grid_id,]) # m
    maxrad <- 3e6
    drast <- exp(log(0.01) * (dist_to_col / maxrad))
    drast[dist_to_col > maxrad] <- 0

    grid_nl[, dens := NA]
    grid_nl[, dens := drast[grid_id]]
    grid_nl[, density_norm := dens / sum(dens * as.numeric(water_area))]

    d <- grid_nl[density_norm > 0, .(code.sra = col$code.sra, site = col$loc_new, grid_id, centre.x, centre.y, density_norm,
                                     area_water = water_area,
                                     source = 'Generated distance to colony')]
    
    rbindlist(lapply(1L:4L, function(q) {
        d[, quarter := q]
        return(d)
    }))

}))

dist <- rbind(dist, dist_missing, fill=T)


## ** Turn normalised densities to numbers of birds
dist[pops_mean, ntot := i.ntot, on = c('code.sra' = 'species_code', 'site')]
dist <- dist[!(code.sra %in% c('PUC', 'DPK'))]

stopifnot(!any(is.na(dist$ntot)))
dist[, density := density_norm * ntot]



## ** Copy whole distribution

## *** format: [target species] = [proxy species]
proxies <- c(DIP = 'DIQ') # Southern royal: use Northern's

for (i in seq_along(proxies)) {
    sp_targ <- names(proxies)[i]
    sp_proxy <- as.character(proxies[i])
    d <- dist[code.sra == sp_proxy]
    d <- d[, .(density = sum(density)), .(code.sra, quarter, grid_id, centre.x, centre.y, area_water)]
    d[, `:=`(code.sra = sp_targ,
             site = 'Proxy distribution',
             source = 'Proxy whole distribution')]
    ## **** Normalise
    d[, density_norm := density / sum(density * as.numeric(area_water))]
    ## **** Merge initial colonies and add ntot
    pm <- pops_mean[species_code == sp_targ,
                    .(site = 'Proxy distribution', ntot = sum(ntot), p_world = 1),
                    .(species_code)]
    pops_mean <- rbind(pops_mean[species_code != sp_targ],
                       pm, fill=T)
    d[pops_mean, ntot := i.ntot, on = c('code.sra' = 'species_code')]
    d[, density := density_norm * ntot]
    dist <- rbind(dist[code.sra != sp_targ], d, fill=T)
}



## ** Aggregate species for each species groups

spcode='ALZ'
group_species <- sapply(sppchar[base_species == F, code.sra], function(spcode) {
    current <- spcode
    allchildren <- NULL
    cont <- T
    while (cont) {
        current <- sppchar[parent.code %in% current, code.sra]
        allchildren <- c(allchildren, current)
        if (!length(current)) {
            cont <- F
        }
    }
    return(sppchar[code.sra %in% allchildren & base_species == T, code.sra])
})

gp='ALZ'
group_dist <- rbindlist(lapply(names(group_species), function(gp) {
    gd <- dist[code.sra %in% group_species[[gp]]]
    gd <- gd[, .(code.sra = gp, site = 'Combined', density = sum(density)), .(centre.x, centre.y, grid_id, quarter)]
}))

dist <- rbind(dist[, .(code.sra, site, quarter, centre.x, centre.y, grid_id, density)], group_dist, fill=T)

dist[sppchar, taxa := i.Taxa, on = 'code.sra']

fwrite(dist, '../generated/gridded-distribution-by-species.csv')


