all: generated/pst-samples.rdata

.first:
	mkdir -p generated  &&  touch .first

## * Demographic parameters

generated/@estimates_table_final.csv: \
		.first \
		../input-data/species-demography/demographic-parameters.csv \
		../input-data/acap-species-list.csv \
		pst_parameter_chosen_values.r \
		prepare_estimates.r
	Rscript prepare_estimates.r

generated/dem-samples.rdata: \
		generated/@estimates_table_final.csv \
		../input-data/species-demography/colony-population-sizes-samples.rds \
		../input-data/acap-species-list.csv \
		create_samples.r
	Rscript create_samples.r

generated/pst-samples.rdata: \
		generated/dem-samples.rdata \
		pst-calc.r
	Rscript pst-calc.r

clean:
	rm -fr generated
