suppressPackageStartupMessages({
    library(data.table)
    library(sf)
    library(ggplot2)
    library(RColorBrewer)
    library(rgdal)
})
options(scipen=100)

source('../../../functions.r')
source('../../../config.r')

load('../../generated/model/model-data-image.rdata', v=F)
load('../../../input-data/grid/grids.rdata')

eff.obs <- overlaps0[type == 'observed' & taxa == taxa[1],
                     .(fishery_group, grid_id, year, quarter, hooks)]

eff.tot <- overlaps0[type == 'total' & taxa == taxa[1],
                     .(fishery_group, grid_id, year, quarter, hooks)]


## ** Important fisheries
fgs <- eff.obs[, unique(fishery_group)]

eff.obs[, hooks := as.numeric(hooks)] # to avoid warning of sum(hooks) being too large for integers
eff.tot[, hooks := as.numeric(hooks)] # to avoid warning of sum(hooks) being too large for integers

eff.obs[, hooks := hooks / 1000]
eff.tot[, hooks := hooks / 1000]

## ** Annual average
eff.tot.yravg <- eff.tot[, .(hooks = sum(hooks) / length(YEARS_PRED)), .(fishery_group, quarter, grid_id)]

## *** All quarters
eff.tot <- rbind(eff.tot.yravg,
                 eff.tot.yravg[, .(quarter = 0, hooks = sum(hooks)), .(fishery_group, grid_id)],
                 fill = T)

## * Gridded distributions

pal <- c("#E3F2FD","#BBDEFB","#90CAF9","#64B5F6","#42A5F5","#2196F3","#1E88E5","#1976D2","#1565C0","#0D47A1")

efftype='tot'
plot_eff <- function(efftype) {
    cat('\n===', efftype, '\n')
    eff <- get(sprintf('eff.%s', efftype))
    fgs <- unique(eff$fishery_group)
    fg=fgs[1]
    for (fg in fgs) {
        cat('\t', fg, '\n')
        quart=0
        for (quart in 0:4) {
            cat(sprintf('\t\tquarter %s\n', ifelse(quart == 0, 'all', quart)))
            griddens <- eff[fishery_group == fg & quarter == quart, .(grid_id, hooks)]
            if (nrow(griddens)) {
                
                g <- map_grid_values(grid_values=griddens,
                                    sprintf("%s\n(x1000 hooks)",
                                            ifelse(efftype == 'tot', 'Mean annual\ntotal effort', 'Total\nobserved effort')),
                                    highlighted.lat=-30, leg.pos = c(0.91, 0.15), pal=pal, colourscale.trans='identity')
                ggsave(sprintf('../assets/map-effort_%s_%s_quarter%i.png',
                               efftype, slugify(fg), quart),
                       width = 7, height = 7, dpi=100)

            }
        }
    }
}

plot_eff('tot')
