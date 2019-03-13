library(shinystan)  # Quick view of Stan output on web browser
load('assets/model-results.rdata')

launch_shinystan(fit) #launch shinyapp for quick view of stan output
