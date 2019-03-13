library(data.table)

load('../../generated/model/model-data-image.rdata', v=F)
spp <- fread('../../acap-species-list.csv')

source('../../../functions.r')
source('../../../labs.r')

spgs <- setdiff(unique(spp$vulgroup), '')

ov <- overlaps[type == 'observed']

cap <- rbind(
    ov[, .(captures = sum(captures)), .(species_group, species_code, quarter)],
    ov[, .(quarter='all', captures = sum(captures)), .(species_group, species_code)])

cap[, species_group := factor(species_group, levels = spgs)]
cap[, sg_name := as.factor(labs[as.character(species_group)])]
cap[, sg_name := reorder(sg_name, as.numeric(species_group))]

cap[spp, sp_name := upper1st(i.cname), on = c('species_code' = 'code.sra')]

capw <- dcast(cap, sg_name + sp_name ~ quarter, value.var = 'captures')
capw <- capw[all > 0]

fwrite(capw, '../assets/table-observed-captures.csv')
