library(data.table)
library(ggplot2)
library(colorspace)

source('../../../functions.r')
source('../../../labs.r')

spplist <- fread('../../acap-species-list.csv')

MAX_RISKRATIO <- 10


## * Summary table

riskratios <- fread('../../generated/risk-ratios-summary.csv')

setorder(riskratios, -rr_p_ov_1, -rr_med)

riskratios <- riskratios[, .(cname = upper1st(cname),
               pst_mean = format(round(pst_mean, 1), big.mark = ' '), pst_ci,
               apf_mean = format(round(apf_mean, 1), big.mark = ' '), apf_ci,
               rr_med = round(rr_med, 2), rr_ci,
               rr_p_ov_1 = round(rr_p_ov_1, 2))]

fwrite(riskratios, '../assets/risk-ratios-summary.csv')


## * Violin plot

risksamples <- readRDS('../../generated/risk-ratios.rds')

get_dens <- function(x) {
    d <- density(x, from = quantile(x, 0.025), to = pmin(quantile(x, 0.975), MAX_RISKRATIO))	
    return(data.table(x = as.numeric(d$x), y = as.numeric(d$y)))	
}

dens <- risksamples[, get_dens(riskratio), .(cname, species_code)]
dens[, y := y / max(y), .(cname)]

dd <- dens[, rbind(.SD, .SD[.N:1, .(x = x, y = -y)], .SD[1]), .(cname, species_code)]

dd[spplist, spgroup := labs[i.vulgroup], on = c('species_code' = 'code.sra')]

stats <- risksamples[, .(rr_med = median(riskratio), p_rr_ov = mean(riskratio > 1)), .(cname, species_code)]

## ** Order taxa
spp <- stats[order(-p_rr_ov, -rr_med), cname]
dd[, lab := factor(as.character(cname), levels = spp)]
    
## ** Median line
dd[stats, med := i.rr_med, on = c('cname')]
med_lims <- dd[, .SD[which.min(abs(med - x))], .(cname, species_code, spgroup)]

outname <- '../assets/plot-risk-ratios.png'

g <- ggplot(dd, aes(x=x, y=y, fill=spgroup)) +
    geom_polygon(alpha = 0.6) +
    geom_segment(data = med_lims, aes(x = med, xend = med, y = y, yend = -y, colour=spgroup)) +
    scale_fill_discrete_qualitative(palette = "Harmonic", name = NULL) +
    scale_colour_discrete_qualitative(palette = "Harmonic", name = NULL) +
    scale_x_continuous(expand = c(0,0), limits = c(0, 1.02*MAX_RISKRATIO)) +
    facet_wrap(~ lab, scales = 'free_y', ncol = 1, strip.position = 'left') +
    labs(x = 'Risk ratio', y = NULL) +
    theme_minimal() +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.text.y =element_text(angle=180, hjust = 1),
          legend.position = c(0.8, 0.2),
          legend.background = element_rect(fill='#EEEEEE', colour=NA),
          legend.margin = margin(0,2,2,2, 'mm'),
          panel.spacing.y=unit(0, 'mm'))

ggsave(outname, g, width = 8, height = 6)

cat('Plot saved as `', outname, '`\n')
