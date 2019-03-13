library(data.table)

source('../../functions.r')

load('../generated/model/mcmc-results.rdata', v=F)

spp <- fread('../acap-species-list.csv')[base_species == T]
setkeyv(spplist, 'code.sra')


c('wasg', 'DIX', 'South Georgia'), 
c('wasg', 'DIX', 'South Georgia'), 

mc_attributes[vartype=='apf_t' & role=='quarter0_spatial2' & species_code=='DIX', .N, .(species_code, site)]
mcsel <- mcmc[variable %in% mc_attributes[role == 'flag_quarter0_spatial0' & vartype == 'apf_t', variable]]
mcsel <- merge(mcsel, mc_attributes, by = 'variable', all.x=T, all.y=F)
apf_samples <- mcsel[, .(apf = sum(value)), .(species_code, sample)]

apf_samples[spplist,
            `:=`(cname = upper1st(i.cname),
                 sname = i.Taxa),
            on = c('species_code' = 'code.sra')]

saveRDS(apf_samples, '../generated/apf-samples.rds')

apf_sp_summ <- apf_samples[
  , .(mean     = mean(apf),
      median   = median(apf),
      sd       = sd(apf),
      lcl      = quantile(apf, 0.025, names=F),
      ucl      = quantile(apf, 0.975, names=F),
      nsamples = .N), .(species_code)]

apf_sp_summ[spplist,
            `:=`(common_name = upper1st(i.cname),
                 sname = i.Taxa),
            on = c('species_code' = 'code.sra')]
setorder(apf_sp_summ, -mean)

fwrite(apf_sp_summ, '../generated/apf-summary.csv')

## * Aggregated spatial APF
#TBC
mcsel <- mcmc[variable %in% mc_attributes[role == 'flag_quarter0_spatial1' & vartype == 'apf_t', variable]]
mcsel <- merge(mcsel, mc_attributes, by = 'variable', all.x=T, all.y=F)
apf_samples <- mcsel[, .(apf = sum(value)), .(grid_id)]


## ** Risk

risk_ratio <- merge(pst_samples, apf_samples, by = c('species_code', 'sample'))
risk_ratio[, riskratio := apf / pst]


## ** Summarise
risk_summ <- risk_ratio[, .(pst_mean = mean(pst),
                            pst_lcl  = quantile(pst, 0.025, names = F),
                            pst_ucl  = quantile(pst, 0.975, names = F),
                            apf_mean = mean(apf),
                            apf_lcl  = quantile(apf, 0.025, names = F),
                            apf_ucl  = quantile(apf, 0.975, names = F),
                            rr_med  = median(riskratio),
                            rr_lcl   = quantile(riskratio, 0.025, names = F),
                            rr_ucl   = quantile(riskratio, 0.975, names = F),
                            rr_p_ov_1 = mean(riskratio > 1)), species_code]
risk_summ[, `:=`(pst_ci = sprintf('%s--%s',
                                  format(round(pst_lcl,1), big.mark=' ', trim=T),
                                  format(round(pst_ucl,1), big.mark=' ', trim=T)),
                 apf_ci = sprintf('%s--%s',
                                  format(round(apf_lcl), big.mark=' ', trim=T),
                                  format(round(apf_ucl), big.mark=' ', trim=T)),
                 rr_ci = sprintf('%s--%s',
                                 format(round(rr_lcl,2), big.mark=' ', trim=T),
                                 format(round(rr_ucl,2), big.mark=' ', trim=T)))]

risk_summ <- merge(spplist[, .(species_code=code.sra, sname=Taxa, cname)], risk_summ, by = 'species_code')


saveRDS(risk_ratio, '../generated/risk-ratios.rds')
fwrite(risk_summ, '../generated/risk-ratios-summary.csv')
    
