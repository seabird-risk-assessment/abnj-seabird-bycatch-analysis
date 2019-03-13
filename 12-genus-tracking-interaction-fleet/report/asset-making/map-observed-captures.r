suppressPackageStartupMessages({
    library(data.table)
})
load('../../generated/model/model-data-image.rdata', v=F)

source('../../../config.r')
source('../../../labs.r')

load('../../../input-data/grid/grids.rdata')

spplist  <- fread('../../acap-species-list.csv')
source('../../../functions.r')

ov <- overlaps[type == 'observed']
ov[, sgroup := ifelse(base_species == T, species_group, species_code)]

ov[, sum(captures), sgroup]

captures <- overlaps[type == 'observed' & base_species == T,
                     .(captures = sum(captures)),
                     .(species_group, grid_id)]


pal <- c("#FEEDB0","#F9C58B","#F39E6E","#EA785A","#D85455","#BE375A","#9D2461","#791A60","#531653","#2F0F3E")

spg=captures$species_group[1]
for (spg in unique(captures$species_group)) {
    cat(spg, '\n')
    d <- captures[species_group == spg, .(grid_id, captures)]
    g <- map_grid_values(grid_values=d, legend.title = 'Observed\ncaptures', pal=pal, colourscale.trans='identity')
    ggsave(sprintf('../assets/map-observed-captures_%s.png', spg), g, width = 7, height = 7, dpi=100)
}

