function [fdot,gdot] = fdot_and_gdot(mu,chi,r,r0,alpha)
% This function calculates the time derivatives of the Lagrange f and g
% coefficients.

% mu = gravitational parameter [km3/s2]
% alpha = reciprocal of the semimajor axis [1/km]
% r0 = radial position at time t0 [km]
% t = time elapsed since r0 [s]
% r = radial position after time t [km]
% chi = universal anomaly after time t [km^0.5]
% fdot = Lagrange f coefficient [1/s]
% gdot = Lagrange g coefficient [ND]

z = alpha*chi^2;
[S,C] = stumpff(z);

fdot = sqrt(mu)/r/r0*(z*S-1)*chi;
gdot = 1 - chi^2/r*C;

end