// Nonlinear Ramsey problem: optimal growth with log utility.
var c k a;
varexo e;
parameters alpha delta rho beta;
alpha = 0.33; delta = 0.025; rho = 0.9; beta = 0.99;

model;
  k = exp(a) * k(-1)^alpha + (1 - delta) * k(-1) - c;
  a = rho * a(-1) + e;
end;

initval;
  a = 0; k = 28; c = 2.3;
end;

shocks;
  var e; stderr 0.01;
end;

planner_objective log(c);
ramsey_model(planner_discount = beta, instruments = (c));
stoch_simul(order = 1, irf = 20);
