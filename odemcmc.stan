functions {
	real[] m1ode(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {        	real dydt[5];        	dydt[1] = theta[1] * y[1] * (1-y[1]-y[2]-y[3]-y[4]-y[5]) - theta[2] * y[1];        	dydt[2] = theta[2] * y[1] + theta[3] * y[2] * (1-y[1]-y[2]-y[3]-y[4]-y[5]) - theta[4] * y[2];
		dydt[3] = theta[4] * y[2] - theta[5] * y[3];
		dydt[4] = theta[6] * y[4] * (1-y[1]-y[2]-y[3]-y[4]-y[5]) - theta[7] * y[4];
		dydt[5] = theta[7] * y[4] - theta[8] * y[5];        	return dydt;
}      
}data {	int<lower=1> T;      	real R[T];      	real t0;      	real ts[T];}transformed data {      	real x_r[0];      	int x_i[0];
}parameters {      	real<lower=0> sigma;
      	real<lower=0,upper=0.2> y0[5];      	real theta[8];}
transformed parameters {
	real<lower=0,upper=1> z;
	z = y0[1] + y0[2] + y0[3] + y0[4] + y0[5];
}
model {      	real y_hat[T,5];
      	real R_hat[T];      	sigma ~ uniform(0,1);
      	theta ~ normal(0.5,1);
	y0 ~ uniform(0,0.2);      	y_hat = integrate_ode_rk45(m1ode, y0, t0, ts, theta, x_r, x_i);
	
      	for (t in 1:T) {
		R_hat[t] = y_hat[t,5]/(y_hat[t,5] + 2 * y_hat[t,3]);        	//R[t] ~ lognormal(log(R_hat[t]) - (sigma^2)/2 , sigma);
		R[t] ~ normal(R_hat[t], sigma);
      	}}
