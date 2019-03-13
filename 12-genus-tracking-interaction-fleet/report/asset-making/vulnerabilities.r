## setwd('stan')
library(data.table)
library(ggplot2)

source('../../functions.r')
source('../../data/parameters.r')

spplist  <- fread('../../data/species-info/species-list.csv')
setkey(spplist, code)
nspp <- nrow(spplist)

modeldat <- new.env()
load('../../generated/data.rdata', modeldat)

load('assets/mcmc-results.rdata')

plot_prior_posterior <- function(mc, mainlab=NULL, sublab=NULL, panel_scales = 'free_y', xlims = NULL, trans_x = NULL) {
    g <- ggplot(mc) +
        geom_density(aes(x = value, y = ..scaled.., fill = type), size = 0.1, colour = 'black') +
        facet_wrap(~ label, scales = panel_scales) +
        scale_x_continuous(limits = xlims, trans = trans_x) +
        scale_fill_manual(name = NULL, values = c(Prior = '#55555555', Posterior = "#238B4599")) +
        labs(x = 'log(Vulnerability)', y = 'Density') +
        theme_minimal() +
        theme(panel.grid.major.y=element_blank(),
              panel.grid.minor.y=element_blank(),
              strip.text=element_text(size = 6),
              axis.text = element_text(size = 6))
    if (!is.null(mainlab)) {
        if (!is.null(sublab)) {
            g <- g + ggtitle(mainlab, sublab)
        } else {
            g <- g + ggtitle(mainlab)
        }
    }
    return(g)
}


## * Table of estimates

qs <- mcsumm[grepl('^q', variable, perl=T) & !grepl('prior', variable, perl=T)]

qs[grepl('^q0', variable, perl=T),
   par_lab := 'Intercept']
qs[grepl('^q_f0', variable, perl=T),
   par_lab := 'Fishery group']
qs[grepl('^q_g0', variable, perl=T),
   par_lab := 'Species group']
qs[grepl('^q_gf0', variable, perl=T),
   par_lab := 'Species x fishery group']

qs[, fg_lab := flags[as.character(fishery_group)]]

qs[species_group == 'ANT', sp_lab := 'Antipodeans']
qs[species_group == 'RYL', sp_lab := 'Royals']
qs[species_group == 'WND', sp_lab := 'Wanderers']
qs[species_group == 'DIP', sp_lab := 'S. royal']
qs[species_group == 'DIQ', sp_lab := 'N. royal']
qs[species_group == 'DIW', sp_lab := 'Gibson\'s']
qs[species_group == 'DQS', sp_lab := 'Antipodean']
qs[species_group == 'TAM', sp_lab := 'Tristan + Amst.']

qs[vartype %in% 'q_gf0', sp_lab := upper1st(spplist[match(species, code), common_name])]
qs[vartype %in% 'q_gf0', lab := sprintf('%s in %s', sp_lab, fg_lab)]
qs[vartype %in% 'q_g0', lab := sp_lab]
qs[vartype %in% 'q_f0', lab := fg_lab]

qs_export <- qs[, .(par_lab, lab, mean, lcl, ucl, min, max)]
qs_export[, par_lab := factor(par_lab, levels = c('Intercept', 'Fishery group', 'Species group', 'Species x fishery group'))]
setorder(qs_export, par_lab, -mean)


fwrite(qs, 'assets/vulnerabilities-all.csv')
fwrite(qs_export, 'assets/vulnerabilities.csv')



## * Plot of prior/posterior

## ** Posteriors
mc <- mcmc[grepl('^q', variable, perl=T) & !grepl('prior', variable, perl=T)]
mc <- merge(mc, mc_attributes, by = c('variable'), all.x=T, all.y=F)
mc[, type := 'Posterior']

## ** Priors
priors <- mcmc[grepl('^q', variable, perl=T) & grepl('prior', variable, perl=T)]
priors_d <- rbindlist(lapply(qs$variable, function(v) {
    if (grepl('^q0|^q_f0|^q_g0', v)) {
        p <- priors[variable == 'q_prior']
    } else {
        p <- priors[variable == 'q_gf_prior']
    }
    p[, variable := v]
    p[, type := 'Prior']
    return(p)
}))
priors_d <- merge(priors_d, mc_attributes, by = c('variable'), all.x=T, all.y=F)


## ** Labels
prior_posterior <- rbind(priors_d, mc)

prior_posterior[grepl('^q0', variable, perl=T),
   par_lab := 'Intercept']
prior_posterior[grepl('^q_f0', variable, perl=T),
   par_lab := 'Fishery']
prior_posterior[grepl('^q_g0', variable, perl=T),
   par_lab := 'Species']
prior_posterior[grepl('^q_gf0', variable, perl=T),
   par_lab := 'Species x fishery']

prior_posterior[, fg_lab := flags[as.character(fishery_group)]]

prior_posterior[species_group == 'ANT', sp_lab := 'Antipodeans']
prior_posterior[species_group == 'RYL', sp_lab := 'Royals']
prior_posterior[species_group == 'WND', sp_lab := 'Wanderers']
prior_posterior[species_group == 'DIP', sp_lab := 'S. royal']
prior_posterior[species_group == 'DIQ', sp_lab := 'N. royal']
prior_posterior[species_group == 'DIW', sp_lab := 'Gibson\'s']
prior_posterior[species_group == 'DQS', sp_lab := 'Antipodean']
prior_posterior[species_group == 'TAM', sp_lab := 'Tristan + Amst.']

prior_posterior[, lab := paste(na.omit(c(sp_lab, fg_lab)), collapse=' in '), by = 1:nrow(prior_posterior)]
prior_posterior[, label := paste(setdiff(c(par_lab, lab), ''), collapse='\n'), by = 1:nrow(prior_posterior)]
prior_posterior[, label := factor(label, levels = unique(label))]

prior_posterior[, type := factor(type, levels = c('Prior', 'Posterior'))]

g <- plot_prior_posterior(prior_posterior, panel_scales='free', trans_x = 'log')

ggsave('assets/vulnerabilities_comparison-prior-posterior.png', g, width = 9, height = 7)

