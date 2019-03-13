library(data.table)

eff <- fread('../../../sra-southern-hemisphere/data/nz-captures-and-effort/nz-effort0.csv')
oeff <- eff[type == 'observed']


cap <- fread('../../../sra-southern-hemisphere/data/nz-captures-and-effort/nz-captures0.csv')

