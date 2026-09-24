// Linear model with two-period leads and lags and a lagged shock.
var x y w;
varexo e1 e2;
parameters a1 a2 b1 b2 c1;
a1 = 0.5; a2 = 0.2; b1 = 0.3; b2 = 0.1; c1 = 0.4;

model(linear);
  x = a1 * x(-1) + a2 * x(-2) + e1 + 0.5 * e1(-1);
  y = b1 * y(+1) + b2 * y(+2) + x + w;
  w = c1 * w(-1) + 0.3 * x(-1) + e2;
end;

shocks;
  var e1; stderr 1;
  var e2; stderr 0.5;
end;

stoch_simul(order = 1, irf = 20);
