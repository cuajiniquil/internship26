functions {
  real int_invN(real a, real b, real phi, real r) {
    if (abs(r) < 1e-12) {
      return (b - a) / phi;
    } else {
      return (exp(r * b) - exp(r * a)) / (phi * r);
    }
  }
}

data {
  int<lower=1> n_events;
  vector[n_events] evt_time;
  array[n_events] int<lower=0, upper=1> evt_type; // 1 = sample, 0 = coalescent
}

parameters {
  real<lower=0> phi;
  real r;
}

model {
  // Priors
  phi ~ lognormal(2, 1);
  r   ~ normal(0, 0.5);

  int n_active = 0;
  for (i in 1:n_events) {
    if (evt_type[i] == 1) n_active += 1;
    else n_active -= 1;

    if (i < n_events && n_active >= 2) {
      real k = n_active;
      real kc2 = k * (k - 1) / 2.0;
      real a = evt_time[i];
      real b = evt_time[i + 1];

      target += -kc2 * int_invN(a, b, phi, r);

      if (evt_type[i + 1] == 0) {
        target += log(kc2) - log(phi) + r * b;
      }
    }
  }
}