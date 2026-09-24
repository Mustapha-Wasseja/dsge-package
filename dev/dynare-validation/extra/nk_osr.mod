// Optimal simple rule: Taylor-rule coefficients minimising a loss.
var y pi r u g;
varexo eu eg;
parameters beta sigma kappa rho_u rho_g phi_pi phi_y;
beta = 0.99; sigma = 1; kappa = 0.1; rho_u = 0.5; rho_g = 0.8;
phi_pi = 1.5; phi_y = 0.5;

model(linear);
  y = y(+1) - 1/sigma * (r - pi(+1)) + g;
  pi = beta * pi(+1) + kappa * y + u;
  r = phi_pi * pi + phi_y * y;
  u = rho_u * u(-1) + eu;
  g = rho_g * g(-1) + eg;
end;

shocks;
  var eu; stderr 0.01;
  var eg; stderr 0.01;
end;

osr_params phi_pi phi_y;
osr_params_bounds;
  phi_pi, 1.01, 5;
  phi_y, 0, 2;
end;
optim_weights;
  pi 1;
  y 0.25;
  r 0.1;
end;
osr(opt_algo = 9);
