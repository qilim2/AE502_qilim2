function [r_vec,v_vec] = r0v02rv(r0_vec,v0_vec,deltat,mu)
% This function computes the new state vector [r_vec, v_vec] from the
% initial state vector [r0_vec, v0_vec] and the elapsed time, deltat.

% mu = gravitational paramter in [km3/s2]
% r0_vec = initial position vector in [km]
% v0_vec = initial velocity vector in [km/s]
% deltat = elapsed time in [s]
% r_vec = final position vector in [km]
% v_vec = final velocity vector in [km/s]

% Magnitudes of initial state vectors
r0 = norm(r0_vec);
v0 = norm(v0_vec);

% Initial radial velocity
vr0 = dot(r0_vec,v0_vec)/r0;

% Reciprocal of the semimajor axis using the energy equation
alpha = 2/r0 - v0^2/mu;

% Compute the universal anomaly
chi = UniversalKepler(mu,deltat,r0,vr0,alpha);

% Compute the f and g functions
[f,g] = f_and_g(mu,chi,deltat,r0,alpha);

% Compute final position vector
r_vec = f*r0_vec + g*v0_vec;
r = norm(r_vec);

% Compute derivatives of f and g
[fdot,gdot]= fdot_and_gdot(mu,chi,r,r0,alpha);

% Compute final velocity vector
v_vec = fdot*r0_vec + gdot*v0_vec;


end