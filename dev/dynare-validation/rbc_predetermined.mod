// RBC with endogenous labour, capital declared predetermined
// (k is beginning-of-period capital), steady state from initval only.
var y c k l z;
varexo e;
parameters alpha beta delta psi rho;
alpha = 0.36; beta = 0.99; delta = 0.025; psi = 1.8; rho = 0.9;
predetermined_variables k;

model;
  1/c = beta / c(+1) * (alpha * y(+1) / k(+1) + 1 - delta);
  psi * c / (1 - l) = (1 - alpha) * y / l;
  y = exp(z) * k^alpha * l^(1 - alpha);
  k(+1) = y - c + (1 - delta) * k;
  z = rho * z(-1) + e;
end;

initval;
  z = 0; l = 0.33; k = 11; y = 1.1; c = 0.8;
end;

shocks;
  var e; stderr 0.007;
end;

stoch_simul(order = 1, irf = 20);
