library(data.table)
library(ggplot2)
library(viridis)

options(scipen = 100)

sppdem <- fread('../../data/species-info/demographic-parameters-processed.csv')

load('../../generated/modeldata.rdata')

spplist <- fread('../../data//species-info/species-list.csv')

make_heatmap <- function(dat, xvar, yvar, fillvar, textvar=NULL, textsubvar=NULL,
                 xlab = '', ylab = '', filllab = '', axis.x.lab.rot = 45, axis.x.lab.hadj = 1,
                 nudge_y_main = ifelse(is.null(textsubvar), 0, 0.1), nudge_y_sub = -0.1, size_main = 3, size_sub = 2,
                 trans = 'log2') {
    dat <- copy(dat)
    dat[, xvar := get(xvar)]
    dat[, yvar := get(yvar)]
    dat[, fillvar := get(fillvar)]
    if (is.null(textvar)) {
        dat[, textvar := round(fillvar, 3)]
    } else dat[, textvar := get(textvar)]
    if (is.null(textsubvar)) {
        dat[, textsubvar := '']
    } else dat[, textsubvar := get(textsubvar)]
    g <- ggplot(dat, aes(x = xvar, y = yvar, fill = fillvar)) +
        geom_tile() +
        geom_text(aes(label = textvar), colour = 'white', nudge_y = nudge_y_main, size = size_main) +
        geom_text(aes(label = textsubvar), colour = 'white', nudge_y = nudge_y_sub, size = size_sub) +
        scale_fill_gradientn(name = filllab, colours = c("#B9B9B9", "#E31A1C"), trans = trans) +
        coord_equal(expand=F) +
        theme_bw() +
        theme(panel.grid = element_blank(),
              legend.position='none',
              axis.text.x = element_text(angle = axis.x.lab.rot, hjust = axis.x.lab.hadj),
              plot.margin = unit(rep(0.2, 4), 'cm')) +
        labs(x = xlab, y = ylab)
    return(g)
}



## * Estimated vulnerabilities

load('assets/mcmc-results.rdata', v=T)
vulsumm <- fread('assets/vulnerabilities-all.csv')

vul <- merge(mcmc, vulsumm, by='variable', all.x=F, all.y=T)

fgroups <- vulsumm[vartype == 'q_f0', .N, .(fishery_group, fg_lab)]
spg <- vulsumm[vartype == 'q_g0', .N, .(species_group, sp_lab)]
spp <- vulsumm[vartype == 'q_gf0', .N, .(species, sp_lab)]
spp[, spg := spplist[match(species, code), vul_group]]

q_samples <- rbindlist(lapply(spp[, species], function(sp) {
    rbindlist(lapply(fgroups[, fishery_group], function(fg) {
        q0 <- vul[variable == 'q0']
        qf0 <- vul[vartype == 'q_f0' & fishery_group == fg]
        qg0 <- vul[vartype == 'q_g0' & species_group == spp[species == sp, spg]]
        qgf0 <- vul[vartype == 'q_gf0' & species == sp & fishery_group == fg]
        if (length(unique(c(nrow(q0), nrow(qf0), nrow(qg0), nrow(qgf0)))) != 1)
            stop('Pb selecting samples')
        q <- q0$value * qf0$value * qg0$value * qgf0$value
        data.table(sp = sp, fg = fg, q = q)
    }))
}))

q_summ <- q_samples[, .(mean = mean(q),
                       lcl = quantile(q, 0.025, names = F),
                       ucl = quantile(q, 0.975, names = F)), .(sp, fg)]

q_summ[, sp_lab := spp[match(sp, species), 'sp_lab']]
q_summ[, fg_lab := fgroups[match(fg, fishery_group), 'fg_lab']]
q_summ[, mean_rnd := round(mean, 3)]
q_summ[, ci := sprintf('(%0.3f - %0.3f)', lcl, ucl)]

ggsave('assets/vulnerability-total-matrix.png',
       make_heatmap(q_summ, 'sp_lab', 'fg_lab', 'mean', 'mean_rnd', 'ci'),
       width = 10, height = 6, scale = 0.7)

ggsave('assets/vulnerability-total-matrix-lcl.png',
       make_heatmap(q_summ, 'sp_lab', 'fg_lab', 'lcl', 'mean_rnd', 'ci'),
       width = 10, height = 6, scale = 0.7)

save(q_samples, q_summ, file = 'vulnerabilities_total.rdata')



## * Observed captures / overlap


agg_vul[sppdem[parameter=='population_total'], pop := i.mean, on = 'species']
agg_vul[, ratio := captures / (overlap * pop)]
agg_vul[, spgroup := sprintf('%s-%s', species_group, species)]

agg_vul[, sp_lab := spp[match(agg_vul$species, spp$species), 'sp_lab']]
agg_vul[, fg_lab := fgroups[match(agg_vul$fishery_group, fgroups$fishery_group), 'fg_lab']]

agg_vul[, ratio_rnd := round(ratio, 4)]

ggsave('assets/vulnerability-empirical-matrix.png',
       make_heatmap(agg_vul, 'sp_lab', 'fg_lab', 'ratio', 'ratio_rnd', trans='log1p'),
       width = 10, height = 6, scale = 0.7)
