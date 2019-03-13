library(data.table)

source('../../functions.r')

load('../generated/model/mcmc-results.rdata', v=F)

spplist <- fread('../acap-species-list.csv')[base_species == T]
setkeyv(spplist, 'code.sra')
### APF by genus
mcsel_genus <- mcmc[variable %in% mc_attributes[role == 'quarter0_spatial1' & vartype == 'apf_t', variable]]
mcsel_genus <- merge(mcsel_genus, mc_attributes, by = 'variable', all.x=T, all.y=F)
apf_samples_genus <- merge(mcsel_genus, spplist[, .(species_code=code.sra, genus=vulgroup)], by='species_code', all.x=T)
apf_genus <- apf_samples_genus[, .(apf=sum(value)), .(sample, genus)][, .(amean=mean(apf), ali=quantile(apf, 0.025), aui=quantile(apf, 0.975)), genus]
write.csv(apf_genus, row.names=FALSE, file='../generated/apf-genus.csv')

