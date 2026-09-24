// NK model estimated with a measurement error on r; with interest-rate
// smoothing, a model-local variable and three shocks.
var y pi r g u v;
varexo eg eu ev;
parameters beta sigma kappa phi_pi phi_y rho_r rho_g rho_u rho_v;
beta = 0.99; sigma = 1; kappa = 0.1;
phi_pi = 1.5; phi_y = 0.125; rho_r = 0.8;
rho_g = 0.9; rho_u = 0.7; rho_v = 0.5;

model(linear);
  # rr = r - pi(+1);
  [name = 'IS curve']
  y = y(+1) - 1/sigma * rr + g;
  [name = 'Phillips curve']
  pi = beta * pi(+1) + kappa * y + u;
  r = rho_r * r(-1) + (1 - rho_r) * (phi_pi * pi + phi_y * y) + v;
  g = rho_g * g(-1) + eg;
  u = rho_u * u(-1) + eu;
  v = rho_v * v(-1) + ev;
end;

shocks;
  var eg; stderr 0.01;
  var eu = 0.005^2;
  var ev; stderr 0.0025;
end;



shocks;
  var r; stderr 0.002;
end;
varobs y pi r;
estimated_params;
  rho_g, 0.9;
end;
estimation(datafile = 'data.csv', mode_compute = 0, nograph, plot_priors = 0);
