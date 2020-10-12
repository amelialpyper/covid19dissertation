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
    
      
      return {dS_dt, dI_dt, dR_dt, dD_dt}; //define our SIRD model for stan, equation (2)
  }
}
data {
  int<lower=1> n_days; //number of observed days
  real y0[4]; //number of initial conditions 
  real t0; //initial time point
  real ts[n_days]; //time points observed
  int<lower=1> N; //population
  int<lower=0> newdeaths[n_days];
}
transformed data {
  real x_r[0]; //real variables used to evaluate the function, which only depend on fixed data
  int x_i[1] = { N }; //integer values used to evaluate the function, which only depend on fixed data
}
parameters {
  real<lower=0> gamma;
  real<lower=0> beta;
  real<lower=0> sigma;
  real<lower=0> phi_inv;
 //our model paramaters and the inverse of our overdispersion paramater
}
transformed parameters{
  real y[n_days, 4];
  real phi = 1. / phi_inv;
  
  {
    real theta[3];
    theta[1] = beta;
    theta[2] = gamma; 
    theta[3] = sigma;

    y = integrate_ode_rk45(sir, y0, t0, ts, theta, x_r, x_i); //solving the SIRD ode with initial conditions
  }
}

model {
//defining our prior distributions for our parmameters, informative wide priors so as not to skew the 
//posterior distribution, a narrower priors for italy sigma to allow for mixing of chains.
  //priors
  beta ~ lognormal(log(0.25),0.5);
  gamma ~ lognormal(log(0.1428),0.5);
  sigma ~ lognormal(log(0.001),0.9); //UK
  //sigma ~ lognormal(log(0.001),0.45); //Italy, change as appropriate
  
phi_inv ~ normal(0,5); 
//sampling distribution
//col(matrix x, int n) - The n-th column of matrix x. Here the number of dead people
newdeaths ~ neg_binomial_2(col(to_matrix(y), 4), phi); //equation (6)
 
}
generated quantities {
  real R0 = beta / (gamma + sigma); //equation used to calculate R0, equation (43)
  real recovery_time = 1 / gamma; 
  real pred_newdeaths[n_days];

  pred_newdeaths = neg_binomial_2_rng(col(to_matrix(y), 4), phi); // predicted deaths estimate for the posterior predictive check.

}
