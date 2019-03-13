MODELS = stan-species-random

TYPE = 

MODELPATHS = $(addprefix models/,$(MODELS))
MODELREPORTS = $(addsuffix /model-report.html,$(MODELPATHS))
DATE := $(shell date +%Y-%m-%d) 

RUN ?= docker run -i --rm --net=host --user=$$(id -u):$$(id -g) -v$$(pwd):/work -w /work --env RUN= dragonflyscience/seabird-risk-assessment:2017-10-18

all: generated/pst-samples.rdata \
			generated/overlaps-and-captures.rdata \
			generated/model/.modeldata \
			generated/model/model-fitted.rdata \
			generated/model/mcmc-results.rdata \
			generated/risk-ratios.rds \
			report


generated/.first:
	mkdir -p generated  &&  touch generated/.first


## * Demographic parameters and PST

generated/@estimates_table_final.csv: \
		generated/.first \
		../input-data/species-demography/demographic-parameters.csv \
		acap-species-list.csv \
		initial-data-preparation/pst_parameter_chosen_values.r \
		initial-data-preparation/prepare_estimates.r
	cd initial-data-preparation  &&  Rscript prepare_estimates.r

generated/dem-samples.rdata: \
		generated/@estimates_table_final.csv \
		../input-data/species-demography/colony-population-sizes-samples.rds \
		acap-species-list.csv \
		initial-data-preparation/create_samples.r
	cd initial-data-preparation  &&  Rscript create_samples.r

generated/pst-samples.rdata: \
		generated/dem-samples.rdata \
		initial-data-preparation/pst-calc.r
	cd initial-data-preparation  &&  Rscript pst-calc.r

generated/gridded-distribution-by-species.csv: \
		generated/dem-samples.rdata \
		../input-data/species-distributions/tracking-distributions-densities.rds \
		../input-data/species-distributions/birdlife-range-maps-densities.rds \
		initial-data-preparation/prepare-bird-distributions.r
	cd initial-data-preparation  &&  Rscript prepare-bird-distributions.r

generated/overlaps-and-captures.rdata: \
		generated/gridded-distribution-by-species.csv \
		../input-data/grid/grids.rdata \
		../input-data/effort-captures/effort-total-RFMO.csv \
		../input-data/effort-captures/NZL/observed-effort-captures.csv \
		initial-data-preparation/calculate-overlaps.r
	cd initial-data-preparation &&  Rscript calculate-overlaps.r


## * <! MODEL !>

generated/model/model-data-image.rdata: generated/model/.modeldata

generated/model/.modeldata: \
		../config.r \
		../input-data/cryptic-mortality/p_observable.rdata \
		acap-species-list.csv \
		generated/overlaps-and-captures.rdata \
		model/prepare-model-data.r
	cd model  &&  Rscript prepare-model-data.r  &&  touch ../generated/model/.modeldata

generated/model/model-fitted.rdata: \
		generated/model/.modeldata \
		model/run-model.r \
		model/model.stan
	cd model  &&  Rscript run-model.r

generated/model/mcmc-results.rdata: \
		generated/model/model-fitted.rdata \
		model/get-results.r
	cd model  &&  Rscript get-results.r

generated/risk-ratios.rds: \
		generated/model/mcmc-results.rdata \
		generated/pst-samples.rdata \
		model/risk-calc.r
	cd model  &&  Rscript risk-calc.r


## * <! REPORTING !>
.PHONY: report
report: generated/risk-ratios.rds
	$(MAKE) -C report

### <! Docker image !> ###
.docker: Dockerfile
	docker build -t dragonflyscience/seabird-risk-assessment .  && \
	docker tag dragonflyscience/seabird-risk-assessment dragonflyscience/seabird-risk-assessment:$(DATE) && \
	docker tag dragonflyscience/seabird-risk-assessment dragonflyscience/seabird-risk-assessment:latest && \
	touch .docker

docker-push: .docker
	docker push dragonflyscience/seabird-risk-assessment:$(DATE) && \
	docker push dragonflyscience/seabird-risk-assessment:latest


## * <! Utils !>

clean:
	rm -fr generated

test:
	@echo $(MODELPATHS)
	@echo $(MODELS)
	@echo $(RUN)
