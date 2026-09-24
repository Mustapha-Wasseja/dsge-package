// STEADY_STATE(), a shock lead, a lagged shock and a varexo_det.
var y z c;
varexo e;
varexo_det d;
parameters a rho;
a = 2; rho = 0.5;

model;
  z = rho * z(-1) + e;
  log(y) = log(a) + z * STEADY_STATE(y) + 0.3 * e(+1) + d;
  c = 0.5 * c(+1) + y - STEADY_STATE(y) + e(-1);
end;

steady_state_model;
  z = 0; y = a; c = 0;
end;

shocks;
  var e; stderr 0.1;
end;

stoch_simul(order = 1, irf = 20);
