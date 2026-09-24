// Linear NK model with a zero lower bound handled by OccBin.
var y pi i inot g;
varexo eg;
parameters beta sigma kappa phi_pi phi_y rho_g ilb;
beta = 0.99; sigma = 1; kappa = 0.1; phi_pi = 1.5; phi_y = 0.125;
rho_g = 0.8; ilb = -0.01;

model(linear);
  y = y(+1) - 1/sigma * (i - pi(+1)) + g;
  pi = beta * pi(+1) + kappa * y;
  inot = phi_pi * pi + phi_y * y;
  [name = 'policy', relax = 'zlb']
  i = inot;
  [name = 'policy', bind = 'zlb']
  i = ilb;
  g = rho_g * g(-1) + eg;
end;

occbin_constraints;
  name 'zlb'; bind inot <= ilb; relax inot > ilb;
end;

shocks(surprise);
  var eg; periods 1; values -0.06;
end;

occbin_setup;
occbin_solver(simul_periods = 30);
