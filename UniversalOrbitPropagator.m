function [r_vec,v_vec] = UniversalOrbitPropagator(r0_vec,v0_vec,DeltaT,mu)
% Two-body orbit propagator
% r0_vec = initial position vector
% v0_vec = initial velocity vector
% DeltaT = time elapsed
% mu = gravitational parameter

% This function requires UniversalKepler.m and stumpff.m to run.

% Ensure that r0_vec and v0_vec are column vectors
[r0_rows,r0_cols] = size(r0_vec);
[v0_rows,v0_cols] = size(v0_vec);
if r0_rows < r0_cols
    r0_vec = transpose(r0_vec);
end
if v0_rows < v0_cols
    v0_vec = transpose(v0_vec);
end


r0 = norm(r0_vec); v0 = norm(v0_vec);
%fprintf('r0 = %f, v0 = %f\n',r0,v0);
vr0 = dot(r0_vec,v0_vec)/r0; % Radial velocity
%fprintf('vr0 = %f\n',vr0);
alpha = 2/r0 - (v0^2)/mu; % Reciprocal of SMA
%fprintf('alpha (1/SMA) = %f\n',alpha);
%a = 1/alpha;

chi = UniversalKepler(mu, DeltaT, r0, vr0, alpha);
%fprintf('chi = %f, chi^2 = %f\n', chi, chi^2);
z = alpha*chi^2;
[S,C] = stumpff(z);
%fprintf('S = %f, C = %f\n',S,C);
f = 1 - ((chi^2)/r0)*C;
g = DeltaT - (1/sqrt(mu))*(chi^3)*S;
%fprintf('f = %f, g = %f\n',f,g);

r_vec = f*r0_vec + g*v0_vec;
r = norm(r_vec);

fdot = sqrt(mu)/(r*r0)*(z*S-1)*chi;
gdot = 1-chi^2/r*C;

v_vec = fdot*r0_vec + gdot*v0_vec;
%v = norm(v_vec);

end
