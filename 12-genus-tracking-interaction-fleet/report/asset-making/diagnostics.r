library(data.table)
library(ggplot2)
source('../../../labs.r')

load('../../generated/model/mcmc-results.rdata')

mcsel <- mcmc[variable %in% mc_attributes[vartype %in% c('q0', 'q_f', 'q_g', 'q_f0', 'q_g0', 'q_gf', 'q_gf0'), variable]]
mcsel <- merge(mcsel, mc_attributes, by = 'variable', all.x = T, all.y = F)

mcsel[vartype %in% c('q_g', 'q_g0') , lab := sprintf('Species group - %s', as.character(species_group))]
mcsel[vartype %in% c('q_f', 'q_f0') , lab := sprintf('Fishery group - %s', as.character(fishery_group))]
mcsel[vartype %in% c('q_gf', 'q_gf0') , lab := sprintf('Sp. group x F. group - %s x %s',
                                                       as.character(species_group),
                                                       as.character(fishery_group))]
mcsel[vartype == 'q0'  , lab := 'Intercept']

mcsel[, lab := factor(lab, levels = c('Intercept',
                                      sort(grep('Fishery', unique(lab), val=T)),
                                      sort(grep('Species', unique(lab), val=T)),
                                      sort(grep('Sp\\. group x', unique(lab), val=T))))]

n.labs <- mcsel[, length(unique(lab))]
n.brks.y <- ifelse(n.labs > 30, 2, 5)

brks <- c(1, 3) * rep(10^(0:10), each=2)

g <- ggplot(mcsel, aes(x = ch_sample, y = value, colour = factor(chain), group = factor(chain))) +
    geom_line(size = 0.2) +
    facet_wrap(~ lab, ncol=1, strip.position = 'left', scales = 'free_y') +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(breaks = brks, trans = 'log1p') +
    ## scale_y_continuous(breaks = function(lims) pretty(lims, n=n.brks.y, min.n=2), trans = 'log1p') +
    labs(x = 'MCMC sample', y = NULL) +
    theme_minimal() +
    theme(legend.position='none',
          strip.text.y=element_text(angle=180, hjust = 1, size = 8),
          strip.placement='outside',
          panel.background=element_rect(fill = '#F4F4F4', color = NA),
          panel.grid.minor=element_blank(),
          panel.grid.major.y=element_line(color = 'white'),
          panel.grid.major.x=element_blank(),
          axis.text = element_text(size = 5),
          axis.title.x = element_text(size = 8))

ggsave('../assets/traces_vulnerabilities.png', width = 8, height = 9)

