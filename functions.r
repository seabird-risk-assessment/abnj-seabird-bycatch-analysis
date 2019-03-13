
## * Utilities

upper1st <- function(x) {
    if (is.factor(x))
        x0 <- levels(x) else x0 <- as.character(x)
    nc <- nchar(x0)
    x0[nc > 1] <- sprintf('%s%s', toupper(substr(x0[nc > 1], 1, 1)), substr(x0[nc > 1], 2, nchar(x0[nc > 1])))
    x0[nc == 1] <- toupper(x0[nc == 1])
    if (is.factor(x))
        levels(x) <- x0 else x <- x0
    return(x)
}

slugify <- function(x) {
    gsub(' +', '_', gsub('\'', '', tolower(x)))
}

## * Statistical

estBetaParams <- function(mu, var) {
    ## Estimate alpha and beta of beta distribution from mean and variance
    ## from http://stats.stackexchange.com/questions/12232/calculating-the-parameters-of-a-beta-distribution-using-the-mean-and-variance
    alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
    beta <- alpha * (1 / mu - 1)
    return(params = list(alpha = alpha, beta = beta))
}


## * For PST

## ** Rmax formula from Niel & Lebreton 2008
funNL<- function(x,a,s) (exp((a+s/(x-s))^-1)-x)^2
lmax_nl <- function(a,s) return(optimise(funNL, c(1,2), a=a, s=s, tol=1e-10)$minimum)
Lmax_nl <- function(a,s,usemc=F) {
    if (usemc) {
        library(parallel)
        return(mcmapply(lmax_nl, a, s))
    } else return(mapply(lmax_nl, a, s))
}
Rmax_NL <- function(s,a,...) return(Lmax_nl(a,s,...)-1)


## * For PST

## formula from Neil & Lebreton 2008
funNL<- function(x,a,s) (exp((a+s/(x-s))^-1)-x)^2
lmax_nl <- function(a,s) return(optimise(funNL, c(1,2), a=a, s=s, tol=1e-10)$minimum)
Lmax_nl <- function(a,s,usemc=F) {
    if (usemc) {
        library(parallel)
        return(mcmapply(lmax_nl, a, s))
    } else return(mapply(lmax_nl, a, s))
}
Rmax_NL <- function(s,a,...) return(Lmax_nl(a,s,...)-1)


 ## Survival variation
surv_var_norm_mean_se <- function(n=1000, smean, sse)
    {
    se_beta <- sse/(smean*(1-smean)) ## calculate se(beta(S)) from se(S), uses delta method
    logit_s <- rnorm(n=n, mean=log(smean/(1-smean)), sd=se_beta) ## apply normal variation to logit(s)
    surv <- exp(logit_s)/(1+exp(logit_s)) ## back-transform
    return(surv)
    }
## surv_var_norm_mean_se1 <- Vectorize(surv_var_norm_mean_se0, SIMPLIFY=F)
## surv_var_norm_mean_se <- function(n=1000, smean, sse)
##   do.call('rbind', surv_var_norm_mean_se1(n=1000, smean, sse))

surv_ci0 <- function(smean, sse)
    {
    se_beta <- sse/(smean*(1-smean)) ## calculate se(beta(S)) from se(S)
    logit_s <- log(smean/(1-smean))
    ll <- logit_s-1.96*se_beta
    ul <- logit_s+1.96*se_beta
    ll <- exp(ll) / (1 + exp(ll))
    ul <- exp(ul) / (1 + exp(ul))
    return(data.frame(ll=ll, ul=ul))
    }
surv_ci1 <- Vectorize(surv_ci0, SIMPLIFY=F)
surv_ci <- function(smean, sse)  do.call('rbind', surv_ci1(smean, sse))

surv_meansd_from_ci0 <- function(lcl, ucl)
    {
    logitlcl <- log(lcl/(1-lcl))
    logitucl <- log(ucl/(1-ucl))
    sdlogit <- (logitucl-logitlcl) / (2*1.96)
    meanlogit <- mean(c(logitlcl,logitucl))
    mean <- exp(meanlogit)/(1+exp(meanlogit))
    sd <- sdlogit*(mean*(1-mean))
    return(data.frame(mean=mean, sd=sd))
    }
surv_meansd_from_ci1 <- Vectorize(surv_meansd_from_ci0, SIMPLIFY=F)
surv_meansd_from_ci <- function(lcl, ucl)  do.call('rbind', surv_meansd_from_ci1(lcl, ucl))


ntot2 <- function(Nadtot, surv, afr, nsims=4000) {
    if (length(Nadtot)==1 | is.null(names(Nadtot))) {
        n <- Nadtot
    } else { ## point estimate or vector of samples
	if (all(names(Nadtot)==c('mean','se'))) { ## mean-se
	    n <- surv_var_norm_mean_se(nsims, smean=Nadtot['mean'], sse=Nadtot['se'])
        } else {
            if (prod(names(Nadtot)==c('min','max'))) { ## min-max
                n <- runif(nsims, min=Nadtot['min'], max=Nadtot['max'])
            }
        }
    }
    if (length(surv)==1 | is.null(names(surv))) {
        s <- surv
    } else {
	if (prod(names(surv)==c('mean','se'))) {
	    s <- surv_var_norm_mean_se(nsims, smean=surv['mean'], sse=surv['se'])
        } else {
            if (prod(names(surv)==c('min','max'))) {
                s <- runif(nsims, min=surv['min'], max=surv['max'])
            }
        }
    }
    if (length(afr)==1 | is.null(names(afr))) {
        a <- afr
    } else {
	if (prod(names(afr)==c('mean','se'))) {
	    a <- rnorm(nsims, mean=afr['mean'], sd=afr['se'])
        } else {
            if (prod(names(afr)==c('min','max'))) {
                a <- runif(nsims, min=afr['min'], max=afr['max'])
            }
        }
    }

    ## maximum age
    m=1+log(0.02)/log(s)

    ## Indiv/adult ratio
    r = s^(1-a)

    ## Total population
    ntot = n*r

    if (min(ntot)<0)
    print('===================================>>>  WARNING!!! Minimum of Ntot is negative!')

    return(list(all=ntot, min=quantile(ntot,0.2), s=s, a=a, r=r, m=m))
    ## for min: take the lower 60 percentile, following Wade
}



## * South Pole-centered map of grid values

map_grid_values <- function(grid_values, legend.title='', colourscale.trans='log1p', zeros_to_na=T,
                            leg.pos = c(0.95, 0.15),
                            pal=c("#E3F2FD","#BBDEFB","#90CAF9","#64B5F6","#42A5F5","#2196F3","#1E88E5","#1976D2","#1565C0","#0D47A1"),
                            bg.col = 'grey95', graticule.col = '#B0B0B0', 
                            highlighted.lat = NA, highlighted.lat.col = 'grey', min0=T, grid = grid_noland, ...) {

    suppressPackageStartupMessages({
        library(data.table)
        library(ggplot2)
        library(colorspace)
        library(sf)
    })
    if (ncol(grid_values) != 2)
        stop('`grid_values` needs to have two columns only: grid ids and values')

    if (any(duplicated(grid_values$grid_id)))
        stop('Duplicated `grid_id` values')

    grid_values <- copy(grid_values)
    setDT(grid_values)
    setnames(grid_values, c('grid_id', 'value'))
    
    values <- grid_values[match(grid$grid_id, grid_id), value]
    if (zeros_to_na) {
        values[values == 0] <- NA
    }
    grid$value <- values

    logtrans <- grepl('log', colourscale.trans)
    if (logtrans) {
            
        mx <- max(pretty(c(0, max(grid$value, na.rm=T))))
        z <- round(mx / (10^(-10:10)))
        z <- (-10:10)[min(which(z > 0 & nchar(z) == 1))]
        brks <- c(1, 3) * rep(10^((z-1):z), each=2)

    } else {
        brks <- pretty(grid$value)        
    }
    if (min0) {
        brks <- unique(c(0, brks))
    }

    bb <- st_bbox(grid)
    world <- suppressMessages(st_crop(oneworld, bb))
    world <- st_segmentize(world, 100000)
    
    g <- ggplot() +
        geom_sf(data=grid, aes(fill = value), colour = graticule.col, size = 0.1) +
        scale_fill_gradientn(colours = pal, name = legend.title,
                             na.value='white', limits = c(0, max(brks)),
                             trans = colourscale.trans, breaks = brks)

    if (!is.na(highlighted.lat)) {
        hlat <- st_as_sf(as(SpatialLines(list(Lines(list(Line(cbind(x = seq(0, 360, 1),
                                                                    y = rep(highlighted.lat, 361)))), ID=1)),
                               proj4string=CRS("+init=epsg:4326 +over")), 'SpatialLines'))
        g <- g + geom_sf(dat = hlat, colour = 'grey')
    }

    g <- g + 
        ## ** world
        geom_sf(data = world, fill = '#AAAAAA', colour = NA, size = 0.1) +
        theme_void() +
        theme(plot.margin = unit(rep(3, 4), 'mm'),
              legend.position=leg.pos,
              legend.margin = margin(5, 2, 2, 2, unit='mm'),
              legend.background=element_rect(fill = '#FFFFFF44'),
              legend.title.align = 0,
              legend.text.align = 0,
              legend.key.width = unit(0.5, 'cm'),
              plot.background = element_rect(fill = bg.col)) +
        coord_sf(crs="+proj=ortho +y_0=0 +lon_0=178 +lat_0=-90.0", expand=F, datum=NA)

    return(g)
}



