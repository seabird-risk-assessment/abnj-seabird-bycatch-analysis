suppressPackageStartupMessages({
    library(data.table)
    library(colorspace)
})
source('../../../functions.r')
load('../../../input-data/grid/grids.rdata')

dens <- fread('../../generated/gridded-distribution-by-species.csv')

spp <- fread('../../acap-species-list.csv')[base_species == T]

pal <- sequential_hcl(10, 'Purples 3', rev=T)

sp='DIM'
for (sp in spp$code.sra) {
    cat(sp, '\n')
    d <- dens[code.sra == sp, .(value = sum(density)), .(grid_id)]
    g <- map_grid_values(grid_values=d, legend.title = 'Birds / km^2', pal=pal, colourscale.trans='identity')
    ggsave(sprintf('../assets/distribution-%s.png', sp), g, width = 7, height = 7, dpi=100)

}
