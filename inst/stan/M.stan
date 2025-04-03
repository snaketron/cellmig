data {
  int<lower=0> N;               // number of cells
  int<lower=0> N_well;          // number of wells
  int<lower=0> N_plate;         // number of plates
  int<lower=0> N_plate_group;   // number of plate groups of wells
  int<lower=0> N_group;         // number groups
  vector[N] y;                  // cell velocity for N cells
  int well_id [N];              // well ID
  int plate_id [N_well];        // plate ID
  int plate_group_id [N_well];  // plate group ID
  int group_id [N_plate_group]; // group ID: treatment x dose
  int offset [N_well];          // offset = 1 (use for batch correction)
}

parameters {
  vector [N_plate] alpha_plate;
  vector [N_group] mu_group;
  
  real <lower=0> sigma_bplate;
  real <lower=0> sigma_wplate;
  
  real mu_shape;
  real <lower=0> sigma_shape;
  
  vector [N_well] z_1;
  vector [N_plate_group] z_2;
  vector [N_well] z_3;
}

transformed parameters {
  vector <lower=0> [N_well] shape;
  vector<lower=0> [N_well] mu;
  vector [N_well] mu_well;
  vector [N_plate_group] mu_plate_group;
  
  mu_plate_group = mu_group[group_id] + sigma_bplate * z_2;
  for(w in 1:N_well) {
    if(offset[w]==1) {
      mu_well[w] = alpha_plate[plate_id[w]] + sigma_wplate * z_1[w];
    } else {
      mu_well[w] = alpha_plate[plate_id[w]] + mu_plate_group[plate_group_id[w]] + sigma_wplate * z_1[w];
    }
  }
  mu = exp(mu_well);
  shape = exp(mu_shape + sigma_shape * z_3);
}

model {
  mu_shape ~ normal(2, 1);
  alpha_plate ~ normal(-2, 1);
  mu_group ~ normal(0, 1);
  sigma_bplate ~ normal(0, 0.5);
  sigma_wplate ~ normal(0, 0.5);
  sigma_shape ~ normal(0, 1);
  z_1 ~ std_normal();
  z_2 ~ std_normal();
  z_3 ~ std_normal();

  y ~ gamma(shape[well_id], shape[well_id] ./ mu[well_id]);
}

generated quantities {
  real y_hat_sample [N_well];
  real log_lik [N];

  y_hat_sample = gamma_rng(shape, shape ./ mu);
  for(i in 1:N) {
    log_lik[i] = gamma_lpdf(y[i] | shape[well_id[i]], shape[well_id[i]] ./ mu[well_id[i]]);
  }
}
