data{    
  int N;        
  int y[N];        
  int n[N];        
  vector[N] dose;        
  int prior_only;        
  int debug;        
  vector[2] pri_mu;        
  vector[4] pri_s;        
  real pri_nu;        
}            
parameters{            
  real b0;        
  real bemax;        
  real<lower=0> bed50;        
  real<lower=0> bhill;        
}        
transformed parameters{        
  vector[N-1] tau;        
  vector[N-1] nu;        
  vector[N] eta;        
        
  tau = bhill * ( log(bed50) - log(dose[2:N]) );        
  nu = -log1p_exp(tau);        
        
  if(bhill < 0){        
    eta = append_row(b0 + bemax, b0 + bemax * exp(nu));        
  }        
  if(bhill > 0){        
    eta = append_row(b0, b0 + bemax * exp(nu));        
  }        
  if(debug){        
    print("tau ", tau);        
    print("nu ", nu);        
    print("eta ", inv_logit(eta));        
  }        
} 
model{    
  target += normal_lpdf(b0        | pri_mu[1], pri_s[1]);    
  target += student_t_lpdf(bemax  | pri_nu, 0, pri_s[2]);    
  target += exponential_lpdf(bed50 | pri_s[3]);    
  target += exponential_lpdf(bhill | pri_s[4]);    
    
  if(!prior_only){    
    target += binomial_logit_lpmf(y | n, eta);    
  }    
}
generated quantities{
  vector[N] p;
  int yrep[N];
  for(i in 1:N){
    p[i] = inv_logit(eta[i]);
    yrep[i] = binomial_rng(n[i], p[i]);
  }
}

