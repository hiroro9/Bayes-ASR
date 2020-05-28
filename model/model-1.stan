//model-1
data {
  int<lower=0> N; //the number of sampeling point
  int<lower=0> D; //the number of directions of observed normal strains
  vector[D] en[N]; //vector of normal strains in D directions [micro strain]
  real t[N]; //sampling time [hour]
  matrix[D, D] I; //unit matrix
  vector[3] n[D]; //directional cosine vectors in D directions
}


parameters {
  real<lower=0, upper=100> s11;//stress tensor [MPa]
  real<lower=0, upper=100> s12;
  real<lower=0, upper=100> s13;
  real<lower=0, upper=100> s22;
  real<lower=0, upper=100> s23;
  real<lower=0, upper=100> s33;
  real<lower=0, upper=100> p0; //pore pressure [MPa]
  real<lower=0, upper=10> K;  //bulk modulus [MPa]
  real<lower=0, upper=10> G; //shear modulus [MPa]
  real<lower=0, upper=100> tv; //relaxation time of volumetric defomation [hour]
  real<lower=0, upper=100> ts; //relaxation time of sheare deformation [hour]
  real<lower=0> sigma; //variance of normal strains [micro strain]
}


transformed parameters {
  vector[D] mu[N];
  vector[3] eigen_values;
  matrix[3, 3] eigen_vectors;
  matrix[3, 3] s; 
    s[1, 1] = s11;
    s[1, 2] = s12;
    s[1, 3] = s13;
    s[2, 1] = s12;
    s[2, 2] = s22;
    s[2, 3] = s23;
    s[3, 1] = s13;
    s[3, 2] = s23;
    s[3, 3] = s33;
  for(i in 1:N){
    for(j in 1:D){
      mu[i, j] = (10^6)*((n[j]'*s*n[j] - (1.0/3.0)*trace(s))*(1 - exp(-t[i]/ts))/(2*G*10^4) 
      + ((1.0/3.0)*trace(s) - p0)*(1 - exp(-t[i]/tv))/(3*K*10^4));
    }
  }
  eigen_values = eigenvalues_sym(s);
  eigen_vectors = eigenvectors_sym(s);
}


model {
  for(i in 1:N) {
    en[i] ~ multi_normal(mu[i], sigma*I);
  }
}

