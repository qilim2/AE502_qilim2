function [f,g] = f_and_g(mu, chi, t, r0, alpha)
% This function computes the Lagrange f and g coefficients

% mu = gravitational parameter [km3/s2]
% alpha = reciprocal of the semimajor axis [1/km]
% r0 = radial position at time t0 [km]
% t = time elapsed since r0 [s]
% chi = universal anomaly after time t [km^0.5]
% f = Lagrange f coefficient [ND]
% g = Lagrange g coefficient [s]

z = alpha*chi^2;
[S,C] = stumpff(z);
f = 1 - chi^2/r0*C;
g = t - 1/sqrt(mu)*chi^3*S;

end