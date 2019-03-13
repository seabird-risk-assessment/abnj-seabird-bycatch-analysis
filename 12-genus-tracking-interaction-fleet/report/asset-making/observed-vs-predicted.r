library(ggplot2)
library(xtable)
library(data.table)

source('../../../functions.r')
## source('../../data/parameters.r')
source('../../../labs.r')
## load('../../generated/data.rdata')

load('../../generated/model/model-data-image.rdata')

spplist  <- fread('../../acap-species-list.csv')
setkey(spplist, code.sra)
spplist <- spplist[base_species == T]
nspp <- nrow(spplist)

## sppdem  <- fread('../../data/species-info/demographic-parameters-processed.csv')
## load('../../../generated//effort-and-captures.rdata', v=F)

## capts <- captures

mcsumm <- readRDS('../../generated/model/mcmc-results-summary.rds')
mcsumm[spplist, common_name := i.cname, on = c('species_code' = 'code.sra')]
mcsumm[, fishery := labs[as.character(fishery_group)]]

setorder(mcsumm, vartype, ind1, ind2)



## * Observed captures

cat('\n* Observed captures\n\n')

pred_vs_obs <- cbind(mcsumm[vartype == 'captures_o'],  obs = CAPTURES_O)

pred_vs_obs[, out := (obs < lcl | obs > ucl)]
cat(sprintf('\n\tNumber of strata with observed value outside 95%% c.i. = %i out of %i strata (%0.2f%%)\n',
            pred_vs_obs[, sum(out)], nrow(pred_vs_obs), 100 * pred_vs_obs[, mean(out)]))


brks <- c(0, c(1, 3) * rep(10^(0:2), each=2))
g <- ggplot(pred_vs_obs, aes(x = obs, y = mean, ymin = lcl, ymax = ucl, colour = out)) +
    geom_errorbar() +
    geom_point(size=2) +
    geom_abline(intercept=0, slope=1) +
    scale_x_continuous(trans = 'log1p', breaks=brks) +
    scale_y_continuous(trans = 'log1p', breaks=brks) +
    scale_colour_discrete(name = NULL) +
    theme_minimal() +
    theme(legend.position=c(0.85, 0.32),
          legend.text = element_text(size = 6),
          legend.key.height = unit(0.4, 'cm'),
          legend.background=element_rect(fill=grey(0.95), colour=NA),
          panel.grid.minor=element_blank()) +
    labs(x = 'Observed number of captures', y = 'Predicted number of captures')
ggsave('../assets/observed-vs-predicted-captures-species.png', width = 7, height = 5)



plot_prior_posterior <- function(mc, mainlab=NULL, sublab=NULL) {
    g <- ggplot(mc) +
        geom_density(aes(x = value, y = ..scaled.., fill = type), size = 0.1, colour = 'black') +
        facet_wrap(~ label, scales = 'free_y') +
        scale_x_continuous(limits = c(0, 1)) +
        scale_fill_manual(name = NULL, values = c(Prior = '#55555555', Posterior = "#238B4599")) +
        labs(x = 'Probability', y = 'Density') +
        theme_minimal() +
        theme(panel.grid.major.y=element_blank(),
              panel.grid.minor.y=element_blank())
    if (!is.null(mainlab)) {
        if (!is.null(sublab)) {
            g <- g + ggtitle(mainlab, sublab)
        } else {
            g <- g + ggtitle(mainlab)
        }
    }
    return(g)
}


## * Probabilities

load('../../generated/model/mcmc-results.rdata', v=F)
nsamples <- max(mcmc$sample)

## ** Posteriors
mc <- mcmc[grep('^p_', variable, perl=T)]
mc <- merge(mc, mc_attributes, by = c('variable'), all.x=T, all.y=F)
mc[, lab1 := sprintf('%s%s', ifelse(is.na(fishery_group), '', as.character(fishery_group)),
                     ifelse(is.na(species_group), '', as.character(species_group)))]
mc[, type := 'Posterior']

## ** Re-create priors from initial parameters
priors <- rbindlist(lapply(unique(mc$variable), function(v) {
        data.table(variable = v, sample = 1:1e5, value = rbeta(1e5, 1, 1))
    }))
priors[, type := 'Prior']

priors[unique(mc[, .(variable, lab1, vartype)]), `:=`(lab1 = i.lab1, vartype = i.vartype), on = 'variable']

mc <- rbind(mc[, .(variable, vartype, lab1, type, sample, value)],
            priors[, .(variable, vartype, lab1, type, sample, value)])
mc[, type := factor(type, levels = c('Prior', 'Posterior'))]

mc[vartype=='p_observable', `:=`(label = 'P(Observable)')]

## ** Plots
g <- plot_prior_posterior(mc)
ggsave('../assets/comparison-with-actual_probabilities.png', width = 9, height = 7)

