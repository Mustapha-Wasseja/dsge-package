// Optimal time-consistent policy (discretion) in a linear NK model.
var y pi r u g;
varexo eu eg;
parameters beta sigma kappa lambda rho_u rho_g;
beta = 0.99; sigma = 1; kappa = 0.1; lambda = 0.25; rho_u = 0.5; rho_g = 0.8;

model(linear);
  y = y(+1) - 1/sigma * (r - pi(+1)) + g;
  pi = beta * pi(+1) + kappa * y + u;
  u = rho_u * u(-1) + eu;
  g = rho_g * g(-1) + eg;
end;

shocks;
  var eu; stderr 0.01;
  var eg; stderr 0.01;
end;

planner_objective pi^2 + lambda * y^2;
discretionary_policy(instruments = (r), irf = 20, planner_discount = beta, nograph, noprint);

