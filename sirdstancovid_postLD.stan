functions {
  real[] sir(real t, real[] y, real[] theta, 
             real[] x_r, int[] x_i) {

      real S = y[1];
      real I = y[2];
      real R = y[3];
      real D = y[4];
      
      real N = x_i[1];
      
      real beta = theta[1];
      real gamma = theta[2];
      real sigma = theta[3];
      
      real dS_dt = -beta * I * S / N;
      real dI_dt =  beta * I * S / N - gamma * I - sigma * I;
      real dR_dt =  gamma * I;
      real dD_dt =  sigma * I ;
    
      
      return {dS_dt, dI_dt, dR_dt, dD_dt};
  }
}
data {
  int<lower=1> n_days; //number of observed days
  real y0[4];
  real t0; //initial time point
  real ts[n_days]; //time points observed
  int<lower=1> N; //population
  int<lower=0> newcases[n_days];
}
transformed data {
  real x_r[0];
  int x_i[1] = { N };
}
parameters {
  real<lower=0> gamma;
  real<lower=0> beta;
  real<lower=0> sigma;
  //real<lower=0> R0;
  real<lower=0> phi_inv;
 
}
transformed parameters{
  real y[n_days, 4];
  real phi = 1. / phi_inv;
  
  {
    real theta[3];
    theta[1] = beta;
    theta[2] = gamma; 
    theta[3] = sigma;

    y = integrate_ode_rk45(sir, y0, t0, ts, theta, x_r, x_i);
  }
}

model {
//real lambda[n_days];
  
  //priors
  beta ~ lognormal(log(0.2),0.5);
  gamma ~ lognormal(log(0.07),0.2);
  sigma ~ lognormal(log(0.004),0.1);
  //sigma ~ normal(0, 1);
  phi_inv ~ normal(0,5);
  //R0 ~ normal(1, 2);
  
  // //likelihood
  // for (i in 1:n_days){
  //   lambda[i] = y[i,4]*N;
  //  }
  // //  
  // //sampling distribution
  // //col(matrix x, int n) - The n-th column of matrix x. Here the number of infected people 
  // 
  // newcases ~ poisson(lambda);
  
 newcases ~ neg_binomial(col(to_matrix(y), 4), phi);
 
}
generated quantities {
  //real lambda[n_days];
  real R0 = beta / (gamma + sigma);
  real recovery_time = 1 / gamma;
  real pred_newcases[n_days];
  
  //  for (i in 1:n_days){
  //    lambda[i] = y[i,2]*N;
  //  }
  //  
  // pred_newcases = poisson_rng(lambda);
  pred_newcases = neg_binomial_rng(col(to_matrix(y), 4), phi);

  
}
