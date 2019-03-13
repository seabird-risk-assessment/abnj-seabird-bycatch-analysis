library(data.table)

sppchar <- fread('../acap-species-list.csv')

dat <- fread('../../input-data/species-demography/demographic-parameters.csv')
dat[species.cname == 'New Zealand white-capped albatross', species.cname := 'white-capped albatross']

dat <- dat[!(species.cname %in% c("pink-footed shearwater", "short-tailed albatross", "waved albatross"))]
dat[, code :=  NA_character_]
dat[, code := sppchar[match(tolower(species.cname), tolower(cname)), code.sra]]

stopifnot(dat[, !any(is.na(code))])

source('../../functions.r')

source('pst_parameter_chosen_values.r')


## transform survival parameters and proportion of breeders from % to [0,1]
for (v in c('point.estimate', 'range.min', 'range.max',  'se', 'sd'))
    dat[parameter %in% c('Scurr', 'Sopt', 'PB'), eval(v) := get(v) / 100]

## take standard error if no standard deviation
dat[is.na(sd) & !is.na(se), sd := se]


case_summ <- list()


#####################################################################################
##  Estimate uncertainty for survival and afr if necessary

## *  SURVIVAL
cond <- dat[, parameter=='Scurr' &
             ((is.na(range.min) & !is.na(range.max)) | (!is.na(range.min) & is.na(range.max)))]
if (any(cond))
    stop('ERROR!!! Single bounds for curr survival not treated!')
cond <- dat[, parameter=='Sopt' &
             ((is.na(range.min) & !is.na(range.max)) | (!is.na(range.min) & is.na(range.max)))]
if (any(cond))
    stop('ERROR!!! Single bounds for opt survival not treated!')

## choose sd for point estimate without uncertainty
cond <- dat[, parameter %in% c('Scurr','Sopt') &
             !is.na(point.estimate) & is.na(sd) & is.na(range.min) & is.na(range.max)]
if (any(cond)) {
    case_summ$s_pointest <- cond
    dat[cond, sd := ifelse(data.quality == 'good', ssdgood,
                    ifelse(data.quality == 'medium', ssdmed, ssdpoor))]
}
## estimate + min-max
## get sd from CI...
cond <- dat[, parameter %in% c('Scurr','Sopt') & is.na(sd) & !is.na(range.min) & !is.na(range.max) & 
             (range.type=='CI' | grepl('MARK',range.type))]
if (any(cond)) {
    case_summ$s_eminmax <- cond
    dat[cond, c('point.estimate', 'sd') := surv_meansd_from_ci(range.min, range.max)]
}
## ...and from symmetrical CI (assuming they are real CI)
cond2 <- dat[, parameter %in% c('Scurr','Sopt') & !is.na(point.estimate) & is.na(sd) & !is.na(range.min) &
              !is.na(range.max) & abs((range.min+range.max)/2 - point.estimate)<1e-5 & !cond]
if (any(cond2)) {
    case_summ$s_eci <- cond2
    dat[cond2, c('point.estimate', 'sd') := surv_meansd_from_ci(range.min, range.max)]
}


## * AGE AT FIRST REPRODUCTION
## only min
cond <- dat[, parameter=='A' & !is.na(range.min) & is.na(range.max) & is.na(point.estimate)]
if (any(cond))	{
    case_summ$a_minonly <- cond
    dat[cond, range.max := range.min * eval(a_mult_minonly)]
}
## only max
cond <- dat[, parameter=='A' & is.na(range.min) & !is.na(range.max) & is.na(point.estimate)]
if (any(cond))	{
    case_summ$a_maxonly <- cond
    dat[cond, range.min := range.max * eval(a_mult_maxonly)]
}
## min and point estimate
cond <- dat[, parameter=='A' & !is.na(range.min) & is.na(range.max) & !is.na(point.estimate)]
if (any(cond))	{
    case_summ$a_emin <- cond
    dat[cond, range.max := 2*point.estimate - range.min]
}
## max and point estimate
cond <- dat[, parameter=='A' & is.na(range.min) & !is.na(range.max) & !is.na(point.estimate)]
if (any(cond))	{
    case_summ$a_emax <- cond
    dat[cond, range.min := 2*point.estimate - range.max]
}
## only point estimate
cond <- dat[, parameter=='A' & is.na(range.min) & is.na(range.max) & !is.na(point.estimate)]
if (any(cond)) {
    case_summ$a_pointest <- cond
    dat[cond, c("range.min", "range.max") := list(floor(dat$point.estimate[cond] * eval(a_multmin_meanonly)),
                                                  ceiling(dat$point.estimate[cond] * eval(a_multmax_meanonly)))]
}


## * PROPORTION OF BREEDERS
if (any(dat[, parameter=='PB' & (!is.na(range.min) | !is.na(range.max)),]))
    stop('ERROR!!! min-max for PB non treated!')
## only mean
cond <- dat[, parameter=='PB' & !is.na(point.estimate) & is.na(sd)]
if (any(cond))	{
    case_summ$pb_esd <- cond
    dat[cond, sd := pb_sd_meanonly]
}

case_summ <- as.data.frame(case_summ)
case_summ <- cbind(code=dat$code, parameter=dat$parameter, case_summ)



#####################################################################################
## *  choose S and A

l <- list()
sp=dat$code[1]
for (sp in sort(unique(dat$code))) {
    lst <- list()
    for (param in c('Scurr','Sopt','A')) {
        min <- max <- mean <- sd <- 0
        ds <- dat[code %in% sp & parameter==param]
        if (sum(ds$proxy.sp=='') == nrow(ds))
            use.proxy.sp <- FALSE else use.proxy.sp <- TRUE
        onlyminmax <- prod(!is.na(c(ds$range.min, ds$range.max)))==1
        someminmax <- sum(!is.na(ds$range.min) & !is.na(ds$range.max))>0
        onlymeansd <- prod(!is.na(c(ds$point.estimate, ds$sd)))==1
        somemeansd <- sum(!is.na(ds$point.estimate) & !is.na(ds$sd))>0
        severalmeansd <- sum(!is.na(ds$point.estimate) & !is.na(ds$sd))>1
        ## only min-max
        if (someminmax & !somemeansd) {
            lst[[param]]['min'] <- min(c(ds$range.min, ds$range.max))
            lst[[param]]['max'] <- max(c(ds$range.min, ds$range.max))
        }
        ## only 1 mean/sd
        if (somemeansd & !severalmeansd & !someminmax) {
            cond <- !is.na(ds$point.estimate) & !is.na(ds$sd)
            lst[[param]]['mean'] <- ds$point.estimate[cond]
            lst[[param]]['sd'] <- ifelse(is.na(ds$sd[cond]), ds$se[cond], ds$sd[cond])
        }
        ## mix of min/max and 1 mean/sd
        if (someminmax & somemeansd & !severalmeansd) { # last condition just to make it exclusive from next if 
            cond <- ds[, !is.na(point.estimate) & !is.na(sd)]
            if (param=='A') {
                ds[cond, range.min := point.estimate - 1.96 * sd]
                ds[cond, range.max := point.estimate + 1.96 * sd]
            }
            if (param %in% c('Scurr','Sopt')) {
                ds[cond, c('range.min', 'range.max') := surv_ci(point.estimate, sd)]
            }
            lst[[param]]['min'] <- min(c(ds$range.min, ds$range.max), na.rm=T)
            lst[[param]]['max'] <- max(c(ds$range.min, ds$range.max), na.rm=T)
        }
        ## several mean/sd (can have min/max too)
        if (severalmeansd) {
            cond <- ds[, !is.na(point.estimate) & !is.na(sd)]
            if (param=='A') {
                ds[cond, c('range.min', 'range.max') :=  list(point.estimate - 1.96 * sd,
                                                              point.estimate + 1.96 * sd)]
            }
            if (param %in% c('Scurr','Sopt')) {
                ds[cond, c('range.min', 'range.max') := surv_ci(point.estimate, sd)]
            }
            lst[[param]]['min'] <- min(c(ds$range.min, ds$range.max), na.rm=T)
            lst[[param]]['max'] <- max(c(ds$range.min, ds$range.max), na.rm=T)
        }
        lst[[param]]['proxysp'] <- use.proxy.sp
    }
    l[[sp]] <- lst
}


#####################################################################################
## * proportion breeding

sp=dat$code[1]
for (sp in sort(unique(dat$code))) {
    spc <- dat$code %in% sp
    biennial <- sppchar[code.sra %in% sp, biennial]
    PB <- list(mean=dat[parameter=='PB' & spc, point.estimate], sd=dat[parameter=='PB' & spc, sd])
    if (length(PB$mean) & length(PB$sd)) {
        l[[sp]]$PB <- PB
    } else if (biennial=='Y') {
        l[[sp]]$PB <- PBbi
    } else if (grepl('^N',biennial)) {
        l[[sp]]$PB <- PBan
    } else l[[sp]]$PB <- PBmix
}



#####################################################################################
## * Convert list to final data frame

i=1
dt <- rbindlist(lapply(1:length(l), function(i) {
    sp <- names(l[i])
    data.table(sp        = sp,
               A.mean    = ifelse(length(l[[i]]$A['mean']), l[[i]]$A['mean'], NA),
               A.sd      = ifelse(length(l[[i]]$A['sd']), l[[i]]$A['sd'], NA),
               A.min     = ifelse(length(l[[i]]$A['min']), l[[i]]$A['min'], NA),
               A.max     = ifelse(length(l[[i]]$A['max']), l[[i]]$A['max'], NA),
               A.proxysp = ifelse(length(l[[i]]$A['proxysp']), l[[i]]$A['proxysp'], NA),
               Scurr.mean    = ifelse(length(l[[i]]$Scurr['mean']), l[[i]]$Scurr['mean'], NA),
               Scurr.sd      = ifelse(length(l[[i]]$Scurr['sd']), l[[i]]$Scurr['sd'], NA),
               Scurr.min     = ifelse(length(l[[i]]$Scurr['min']), l[[i]]$Scurr['min'], NA),
               Scurr.max     = ifelse(length(l[[i]]$Scurr['max']), l[[i]]$Scurr['max'], NA),
               Scurr.proxysp = ifelse(length(l[[i]]$Scurr['proxysp']), l[[i]]$Scurr['proxysp'], NA),
               Sopt.mean    = ifelse(length(l[[i]]$Sopt['mean']), l[[i]]$Sopt['mean'], NA),
               Sopt.sd      = ifelse(length(l[[i]]$Sopt['sd']), l[[i]]$Sopt['sd'], NA),
               Sopt.min     = ifelse(length(l[[i]]$Sopt['min']), l[[i]]$Sopt['min'], NA),
               Sopt.max     = ifelse(length(l[[i]]$Sopt['max']), l[[i]]$Sopt['max'], NA),
               Sopt.proxysp = ifelse(length(l[[i]]$Sopt['proxysp']), l[[i]]$Sopt['proxysp'], NA),
               PB.mean   = ifelse(length(l[[i]]$PB$mean), l[[i]]$PB$mean, NA),
               PB.sd     = ifelse(length(l[[i]]$PB$sd), l[[i]]$PB$sd, NA),
               PB.min    = ifelse(length(l[[i]]$PB$min), l[[i]]$PB$min, NA),
               PB.max    = ifelse(length(l[[i]]$PB$max), l[[i]]$PB$max, NA)
               )
}))

## Check completedness
acols <- grep('^A.',names(dt))
scurrcols <- grep('^Scurr.',names(dt))
soptcols <- grep('^Sopt.',names(dt))
pbcols <- grep('^PB.',names(dt))
a <- apply(dt,1,function(x) {
    A <- sum(!is.na(x[acols]))>1
    Scurr <- sum(!is.na(x[scurrcols]))>1
    Sopt <- sum(!is.na(x[soptcols]))>1
    PB <- sum(!is.na(x[pbcols]))>1
    return(c(a=A,scurr=Scurr,sopt=Sopt,pb=PB))
})
a <- as.data.frame(t(a))

b <- apply(a,2,all)
incompl <- names(b[!b])
if (length(incompl)) {
    inc <- which(a[,incompl]==F)
    l[dt$sp[inc]]
}


#############################
## Add common and group names
dt[, cname   := sppchar[match(dt$sp, code.sra), cname]]
dt[, spgroup := sppchar[match(dt$sp, code.sra), Taxa]]

#############################
## Point estimates
for (param in c('A','Scurr','Sopt','PB')) {
    dt[!is.na(get(paste0(param,'.mean'))) & !is.na(get(paste0(param,'.sd'))),
       eval(paste0(param,'.ptest')) := get(paste0(param,'.mean'))]
    dt[!is.na(get(paste0(param,'.min'))) & !is.na(get(paste0(param,'.max'))),
       eval(paste0(param,'.ptest')) := (get(paste0(param,'.min')) + get(paste0(param,'.max')))/2]
    if (any(dt[, is.na(get(paste0(param,'.ptest')))]))
        stop(paste0('ERROR!!! Some NAs for', param))
}


#####################################################################################
##  EXPORT FINAL TABLE
fwrite(dt, '../generated/@estimates_table_final.csv')
