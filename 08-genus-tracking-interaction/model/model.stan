data {
  int<lower=1>  N_SPECIES;
  int<lower=1>  N_FAMILY_GROUP;
  int<lower=1>  N_ROW_O;
  int<lower=1>  N_ROW_T;
  int<lower=1>  N_SPECIES_VUL_GROUP;
  int<lower=1>  SPECIES_O [N_ROW_O];
  int<lower=1>  SPECIES_T [N_ROW_T];
  int<lower=1>  SPECIES_GROUP_O [N_ROW_O];
  int<lower=1>  SPECIES_GROUP_T [N_ROW_T];
  vector<lower=0> [N_ROW_O] OVERLAP_O;
  real<lower=0> OVERLAP_T [N_ROW_T];
  int <lower=0> CAPTURES_O [N_ROW_O];
  
  //additional data//
  int<lower=1>  N_INTERACTION_CAT_O;
  int<lower=1>  N_INTERACTION_CAT_T;
  int<lower=1>  N_FISHERY_GROUP;
  int<lower=1>  N_OBS_FISHERY_GROUP;
  vector<lower=0> [N_ROW_O] Q_SCALE_O;
  vector<lower=0> [N_ROW_T] Q_SCALE_T;
    
  matrix [N_ROW_O, N_SPECIES_VUL_GROUP] DUMMYVAR_SPECIES_GROUP_O;
  matrix [N_ROW_O, N_OBS_FISHERY_GROUP]     DUMMYVAR_FISHERY_GROUP_O;
  matrix [N_ROW_O, N_INTERACTION_CAT_O]   DUMMYVAR_INTERACTION_O;
  int<lower=1>  FISHERY_GROUP_O [N_ROW_O];
  int<lower=1>  FISHERY_GROUP_T [N_ROW_T];
  
  matrix [N_ROW_T, N_SPECIES_VUL_GROUP] DUMMYVAR_SPECIES_GROUP_T;
  matrix [N_ROW_T, N_FISHERY_GROUP]     DUMMYVAR_FISHERY_GROUP_T;
  matrix [N_ROW_T, N_INTERACTION_CAT_T]   DUMMYVAR_INTERACTION_T;

  // Unidentified captures
  int<lower=0> AVES    [N_OBS_FISHERY_GROUP];
  int AVES_N;
  int AVES_FISHERY_O    [AVES_N];
  int AVES_START_O     [AVES_N];
  int AVES_END_O       [AVES_N];
  int AVES_PRIOR       [AVES_N];
  int FAMILY_FISHERY_N;
  int FAMILY_GROUP_O     [N_ROW_O];
  int FAMILY           [FAMILY_FISHERY_N];
  int FAMILY_FISHERY_O [FAMILY_FISHERY_N];
  int FAMILY_FAMILY_O  [FAMILY_FISHERY_N];
  int FAMILY_START_O  [FAMILY_FISHERY_N];
  int FAMILY_END_O [FAMILY_FISHERY_N];
}

parameters{
  vector<lower=0, upper=1> [N_OBS_FISHERY_GROUP]    p_identified_bird;
  matrix<lower=0,upper=1>  [N_OBS_FISHERY_GROUP, N_FAMILY_GROUP] p_identified_family;
  real                                          log_q0;
  real                                          log_q_prior;
  vector [N_SPECIES_VUL_GROUP]                  log_q_g0;
  vector [N_OBS_FISHERY_GROUP]                   log_q_f0;
  real<lower=1e-3>                                           sigma_f;
  real<lower=1e-3>                                          sigma_g;
  real<lower=1e-3>                                          sigma_gf;
  real<lower=1e-3>                                         sigma_prior;
  vector [N_INTERACTION_CAT_O]                    log_q_gf0;
  real                                          log_q_gf_prior;
}


transformed parameters{
  real<lower=0>                                 q0;
  real<lower=0>                                 q_prior;
  vector<lower=0> [N_SPECIES_VUL_GROUP]         q_g0;
  vector<lower=0> [N_OBS_FISHERY_GROUP]             q_f0;
  vector<lower=0> [N_INTERACTION_CAT_O]           q_gf0;
  real<lower=0>                                 q_gf_prior;
  vector          [N_ROW_O]                     log_q_g_o;
  vector          [N_ROW_O]                     log_q_f_o;
  vector          [N_ROW_O]                     log_q_gf_o;
  vector          [N_ROW_O]                     log_q_o;
  vector          [N_ROW_O]                     mu_observable_incidents_o;


  log_q_g_o = DUMMYVAR_SPECIES_GROUP_O * log_q_g0;
  log_q_f_o = DUMMYVAR_FISHERY_GROUP_O * log_q_f0;
  log_q_gf_o = DUMMYVAR_INTERACTION_O * log_q_gf0;

  log_q_o = log_q_g_o  +  log_q_f_o + log_q_gf_o;  //catchability

  mu_observable_incidents_o = Q_SCALE_O  .* (exp(log_q_o) .* OVERLAP_O);
  
  // Back transform for analysis
  q0 =  exp(log_q0);
  q_prior = exp(log_q_prior);
  q_g0 = exp(log_q_g0);
  q_f0 = exp(log_q_f0);
  q_gf0 = exp(log_q_gf0);
  q_gf_prior = exp(log_q_gf0[1]);

}

model{
  
  for (i in 1:AVES_N){
 	p_identified_bird[i] ~ beta(1, AVES_PRIOR[i]);
  }

  for (i in 1:N_OBS_FISHERY_GROUP){
  	for (j in i:N_FAMILY_GROUP){
	     p_identified_family[i, j] ~ beta(1, 1);
	}
  }
  
  // Intercept
  log_q0 ~ normal(0, 1); 
  sigma_f ~ normal(0, 1);
  sigma_g ~ normal(0, 1);
  sigma_gf ~ normal(0, 1);
  sigma_prior ~ normal(0, 1);
  log_q_f0 ~ normal(0, sigma_f);
  log_q_g0 ~ normal(0,  3);
  log_q_gf0 ~ normal(0,  sigma_gf);
  log_q_prior ~ normal(0, sigma_prior);

  for(i in 1:N_ROW_O){
    CAPTURES_O[i] ~ poisson(mu_observable_incidents_o[i] * p_identified_bird[FISHERY_GROUP_O[i]] * p_identified_family[FISHERY_GROUP_O[i], FAMILY_GROUP_O[i]]);
  }
 
  for (f in 1:FAMILY_FISHERY_N){
	FAMILY[f] ~ poisson(sum(p_identified_bird[FAMILY_FISHERY_O[f]] * 
		(1 - p_identified_family[FAMILY_FISHERY_O[f], FAMILY_FAMILY_O[f]]) * 
		mu_observable_incidents_o[FAMILY_START_O[f] : FAMILY_END_O[f]]));
  }

  for (a in 1:AVES_N){
	AVES[a] ~ poisson(sum((1 - p_identified_bird[AVES_FISHERY_O[a]]) * mu_observable_incidents_o[AVES_START_O[a]: AVES_END_O[a]]));
  }
  
}

generated quantities{
  real  q_t [N_ROW_T];
  vector [N_ROW_T] log_q_g_t;
  vector [N_ROW_T] log_q_f_t;
  vector [N_ROW_T] log_q_gf_t;
  vector [N_ROW_O] log_lik;
  int   apf_t[N_ROW_T];
  int   observable_captures_t[N_ROW_T];
  int<lower=1, upper=N_OBS_FISHERY_GROUP>  fishery_rf[N_FISHERY_GROUP];
  int<lower=1, upper=N_OBS_FISHERY_GROUP * N_SPECIES_VUL_GROUP>  interaction_rf[N_FISHERY_GROUP * N_SPECIES_VUL_GROUP];
  vector<lower=0, upper=1> [N_OBS_FISHERY_GROUP]   p_sampling;

  // Observed strata
  int <lower=0> captures_o [N_ROW_O];

  // Generate some random fisheries
  for (f in 1: N_OBS_FISHERY_GROUP){
  	p_sampling[f] = 1.0/N_OBS_FISHERY_GROUP;
	fishery_rf[f] = f;
  }

  for (f in N_OBS_FISHERY_GROUP + 1 : N_FISHERY_GROUP){
	fishery_rf[f] = categorical_rng(p_sampling);
  }

  for (f in 1: N_FISHERY_GROUP){
  	for (s in 1: N_SPECIES_VUL_GROUP) {
		interaction_rf[(f - 1)*N_SPECIES_VUL_GROUP + s] = (fishery_rf[f] - 1) * N_SPECIES_VUL_GROUP + s;
	}
  }
  
  for(i in 1:N_ROW_O){
    captures_o[i] = poisson_rng(mu_observable_incidents_o[i] * p_identified_bird[FISHERY_GROUP_O[i]] * p_identified_family[FISHERY_GROUP_O[i], FAMILY_GROUP_O[i]]);
    log_lik[i] = poisson_lpmf(CAPTURES_O[i] |  mu_observable_incidents_o[i] * p_identified_bird[FISHERY_GROUP_O[i]] * p_identified_family[FISHERY_GROUP_O[i], FAMILY_GROUP_O[i]]);
  }

  // All strata
  log_q_g_t  = DUMMYVAR_SPECIES_GROUP_T * log_q_g0;
  log_q_f_t  = DUMMYVAR_FISHERY_GROUP_T * log_q_f0[fishery_rf];
  log_q_gf_t = DUMMYVAR_INTERACTION_T   * log_q_gf0[interaction_rf];
  
  for(j in 1:N_ROW_T){
    //log_q_o = log_q_g_o  +  log_q_f_o + log_q_gf_o;  //catchability
    q_t[j] = Q_SCALE_T[j] * exp(log_q_g_t[j] + log_q_f_t[j] + log_q_gf_t[j]) * OVERLAP_T[j];  
    //Total incidents
    apf_t[j] = poisson_rng(q_t[j]);
 
  }

}
