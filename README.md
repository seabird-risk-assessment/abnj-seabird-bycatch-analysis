# ABNJ Seabird bycatch analysis

This is the code used for usng the seabird risk assessment method
to analyse seabird bycatch in southern hemisphere (south of 20°S) surface longline fisheries.
The work was carried out as part of the Areas Beyond National Jurisdiction (ABNJ)
global seabird bycatch (GBSA) project. The analysis was carried out at a workshop in 
Kruger National Park, South Africa, during Fenruary 2019. 

The seabird risk assessment method was applied to observed bycatch data collated
from the countries who particpated in the workshop; to seabird distributions
developed by Birdlife International; and to fishing effort data collated from
tuna Regional Fisheries Management Organisations (tRFMOs).  

These input data sets  were not able to be disseminated, and so the code
in this repository does not work from end to end. However, it provides
an example of the application of the risk assessment method. The model
is fitted using Stan and R, and was run on a Linux computer. We expect that
this repository woul donly be useful to people able to read R code, and with
a familiarity of the seabird risk assessment (see the 
[risk assessment project](https://github.com/seabird-risk-assessment/seabird-risk-assessment) 
for an example of a risk assessment applied to test data). 
The steps used for running the models are indicated in the makefiles. 

## The final models

The two models that were included in the reporting both assumed that the
susceptability of seabirds to capture was the same for all birds within
each genus, an both models included an interaction term between the 
genus and the fishery. All captures were assumed to be of dead birds (and
so survivability of live releases was not considered), and no
cryptic mortalities were included in the models. 

The two models differ only in their assumptions about how to
extropolate to fishing by unobserved fleets. In the first model
( in directory `08-genus-tracking-interaction`, the unobserved
fleets are sampled randomly from the observed fleets). In the
second model ( in directory `12-genus-tracking-interaction-fleet`) all
the unobserved fleets are set to a single fleet, which was assumed to have a 
catchability that was similar to other high seas fisheries. 



## Reuse

This code is made available under the [MIT Licence](LICENSE).

