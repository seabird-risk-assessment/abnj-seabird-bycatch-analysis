all: assets report.html

.PHONY: first
first:
	mkdir -p assets

assets: \
		first \
		assets/observed-vs-predicted-captures-species.png \
		assets/traces_vulnerabilities.png \
		assets/.apf-maps \
		assets/.dist-maps \
		assets/.ocap-maps \
		assets/table-observed-captures.csv \
		assets/.risk \
		assets/.apf-maps2 
		#assets/.eff-maps 


assets/observed-vs-predicted-captures-species.png: \
		../../labs.r \
		../generated/model/model-data-image.rdata \
		../acap-species-list.csv \
		../generated/model/mcmc-results-summary.rds \
		../generated/model/mcmc-results.rdata \
		asset-making/observed-vs-predicted.r
	cd asset-making  &&  Rscript observed-vs-predicted.r

assets/traces_vulnerabilities.png: \
		../../labs.r \
		../generated/model/mcmc-results.rdata \
		asset-making/diagnostics.r
	cd asset-making  &&  Rscript diagnostics.r

assets/.apf-maps: \
		../../functions.r \
		../acap-species-list.csv \
		../generated/model/mcmc-results.rdata \
		../generated/model/mcmc-results-summary.rds \
		asset-making/apf-maps.r
	cd asset-making  &&  Rscript apf-maps.r  &&  touch ../assets/.apf-maps

assets/.apf-maps2: \
		../../functions.r \
		../acap-species-list.csv \
		../generated/model/mcmc-results.rdata \
		../generated/model/mcmc-results-summary.rds \
		asset-making/apf-by-site.r
	cd asset-making  &&  Rscript apf-by-site.r  &&  touch ../assets/.apf-maps2

assets/.dist-maps: \
		../../functions.r \
		../generated/gridded-distribution-by-species.csv \
		../acap-species-list.csv \
		asset-making/bird-distributions.r
	cd asset-making  &&  Rscript bird-distributions.r  &&  touch ../assets/.dist-maps

assets/.eff-maps: \
		../../functions.r \
		../../config.r \
		../generated/model/model-data-image.rdata \
		asset-making/effort-maps.r
	cd asset-making  &&  Rscript effort-maps.r  &&  touch ../assets/.eff-maps

assets/.ocap-maps: \
		../../functions.r \
		../generated/model/model-data-image.rdata \
		asset-making/map-observed-captures.r
	cd asset-making  &&  Rscript map-observed-captures.r  &&  touch ../assets/.ocap-maps

assets/table-observed-captures.csv: \
		../generated/model/model-data-image.rdata \
		../acap-species-list.csv \
		../../labs.r \
		asset-making/table-observed-captures.r
	cd asset-making  &&  Rscript table-observed-captures.r

assets/.risk: \
		../acap-species-list.csv \
		../generated/risk-ratios-summary.csv \
		../generated/risk-ratios.rds \
		asset-making/risk-ratios.r
	cd asset-making  &&  Rscript risk-ratios.r  &&  touch ../assets/.risk



html: report.html
report.html: assets \
		custom.css \
		report.Rmd 
	Rscript -e 'rmarkdown::render("report.Rmd", clean=F)'



cleanall:
	rm -fr assets/*.*
