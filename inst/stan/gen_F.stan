data {
    int<lower=0> N_plate;           // number of plates
    int<lower=0> N_group;           // number groups
    int<lower=0> N_well_reps;       // number groups
    int offset;
    // priors
    real prior_alpha_p_M;         // prior mean of alpha_p
    real prior_alpha_p_SD;        // prior SD of alpha_p
    real prior_sigma_bio_M;       // prior mean of sigma_bio
    real prior_sigma_bio_SD;      // prior SD of sigma_bio
    real prior_sigma_tech_M;      // prior mean of sigma_tech
    real prior_sigma_tech_SD;     // prior SD of sigma_tech
    real prior_kappa_mu_M;        // prior mean of kappa_mu
    real prior_kappa_mu_SD;       // prior SD of kappa_mu
    real prior_kappa_sigma_M;     // prior mean of kappa_sigma
    real prior_kappa_sigma_SD;    // prior SD of kappa_sigma
    real prior_mu_group_M;        // prior mean of mu_group
    real prior_mu_group_SD;       // prior SD of mu_group
}

transformed data {
    int<lower=0> N_well = N_plate * N_group * N_well_reps;
}

parameters {}

transformed parameters {}

model {}

generated quantities {
    real y [N_well];
    vector [N_plate] mu_plate_group [N_group];
    vector<lower=0> [N_well] kappa;
    vector<lower=0> [N_well] mu;
    vector [N_well] mu_well;
    real <lower=0> sigma_bio;
    real <lower=0> sigma_tech;
    vector [N_plate] alpha_p;
    vector [N_group] mu_group;
    real kappa_mu;
    real kappa_sigma;
    int well_id;
    
    for(i in 1:N_plate) {
        alpha_p[i] = normal_rng(prior_alpha_p_M, prior_alpha_p_SD);
    }
    for(i in 1:N_group) {
        mu_group[i] = normal_rng(prior_mu_group_M, prior_mu_group_SD);
    }
    sigma_bio = abs(normal_rng(prior_sigma_bio_M, prior_sigma_bio_SD));
    sigma_tech = abs(normal_rng(prior_sigma_tech_M, prior_sigma_bio_SD));
    kappa_mu = normal_rng(prior_kappa_mu_M, prior_kappa_mu_SD);
    kappa_sigma = abs(normal_rng(prior_kappa_sigma_M, prior_kappa_sigma_SD));
    
    well_id = 1;
    for(g in 1:N_group) {
        for(p in 1:N_plate) {
            
            mu_plate_group[g][p] = normal_rng(mu_group[g], sigma_bio);
            
            for(w in 1:N_well_reps) {
                if(g==offset) {
                    mu_well[well_id] = normal_rng(alpha_p[p], sigma_tech);
                }
                if(g!=offset) {
                    mu_well[well_id] = normal_rng(alpha_p[p] + mu_plate_group[g][p], sigma_tech);
                }
                kappa[well_id] = exp(normal_rng(kappa_mu, kappa_sigma));
                well_id = well_id + 1;
            }
        }
        mu =  exp(mu_well);
    }
    
    y = gamma_rng(kappa, kappa ./ mu);
}
