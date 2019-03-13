suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(colorspace)
    library(RColorBrewer)
    library(sf)
})
source('../../../functions.r')

spplist  <- fread('../../acap-species-list.csv')[base_species==T]
setkey(spplist, code.sra)
nspp <- nrow(spplist)

load('../../generated/model/mcmc-results.rdata')
mcsumm <- readRDS('../../generated/model/mcmc-results-summary.rds')

load('../../../input-data/grid/grids.rdata')


## * APF by species

mcsel <- mcmc[variable %in% mc_attributes[role == 'quarter0_spatial2' & vartype == 'apf_t', variable]]
mcsel <- merge(mcsel, mc_attributes, by = 'variable', all.x=T, all.y=F)

## ** By site, non-gridded
apf_site <- mcsel[, .(value = sum(value)), .(species_code, site, sample)][
  , .(grid_id  = 0,
      mean     = mean(value),
      median   = median(value),
      sd       = sd(value),
      lcl      = quantile(value, 0.025, names=F),
      ucl      = quantile(value, 0.975, names=F),
      nsamples = .N),
    .(code = species_code, site)]

apf_site <- apf_site[spplist[base_species == T], on = c('code' = 'code.sra')][
  , .(common_name = upper1st(cname), code, site, grid_id, mean , sd, lcl, ucl)]

fwrite(apf_site, '../assets/apf-by-site.csv')


### ** By site, gridded
#
#apf_sp <- mcsel[, .(value = sum(value)), .(species_code, site, grid_id, sample)][
#  , .(mean     = mean(value),
#      median   = median(value),
#      sd       = sd(value),
#      lcl      = quantile(value, 0.025, names=F),
#      ucl      = quantile(value, 0.975, names=F),
#      nsamples = .N),
#    .(code = species_code, site, grid_id)]
#
#apf_sp <- apf_sp[spplist[base_species == T], on = c('code' = 'code.sra')][
#  , .(common_name = upper1st(cname), code, site, grid_id, mean , sd, lcl, ucl)]
#
#apf_sp <- apf_sp[!is.na(site)]
#
#
### * Maps
#
#pal <- sequential_hcl(10, 'Reds 3', rev=T)
#
#spp <- spplist[base_species == T]
#sp='DIM'
#for (sp in spp$code.sra) {
#    cat(sp, '\n')
#    a0 <- apf_sp[code == sp]
#    for (site1 in unique(a0$site)) {
#        cat('\t', site1, '\n')
#        a <- a0[site %in% site1, .(grid_id, mean)]
#        if (nrow(a) > 1) {
#            g <- map_grid_values(grid_values=a, legend.title = 'APF', pal=pal)
#            ggsave(sprintf('../assets/apf-map_%s_%s.png', sp, slugify(site1)), width = 7, height = 7, dpi=100)
#        }
#    }
#}
#
