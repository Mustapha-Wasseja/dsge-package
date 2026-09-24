// Basic real business cycle model in standard Dynare timing
// (k is end-of-period capital, so production uses k(-1)).
// Import with: read_dynare(system.file("examples", "rbc.mod", package = "dsge"))

var y c k i a;
varexo e;

parameters alpha beta delta rho sigma_e;

alpha   = 0.33;
beta    = 0.99;
delta   = 0.025;
rho     = 0.95;
sigma_e = 0.01;

model;
  1/c = beta / c(+1) * (alpha * exp(a(+1)) * k^(alpha - 1) + 1 - delta);
  y   = exp(a) * k(-1)^alpha;
  k   = i + (1 - delta) * k(-1);
  y   = c + i;
  a   = rho * a(-1) + e;
end;

steady_state_model;
  a = 0;
  k = ((1/beta - 1 + delta) / alpha)^(1 / (alpha - 1));
  y = k^alpha;
  i = delta * k;
  c = y - i;
end;

shocks;
  var e; stderr sigma_e;
end;

steady;
check;
stoch_simul(order = 1, irf = 40) y c k i;
