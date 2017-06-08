functions {
	real[] m1ode(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
		dydt[3] = theta[4] * y[2] - theta[5] * y[3];
		dydt[4] = theta[6] * y[4] * (1-y[1]-y[2]-y[3]-y[4]-y[5]) - theta[7] * y[4];
		dydt[5] = theta[7] * y[4] - theta[8] * y[5];
}      
}
}
      	real<lower=0,upper=0.2> y0[5];
model {
      	real R_hat[T];
      	theta ~ normal(0.5,1);
	y0 ~ uniform(0,0.2);
	
      	for (t in 1:T) {
		R_hat[t] = y_hat[t,5]/(y_hat[t,5] + 2 * y_hat[t,3]);
		R[t] ~ normal(R_hat[t], sigma);
      	}