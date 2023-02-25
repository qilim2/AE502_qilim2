function [v0_vec,v_vec] = LambertSolver_UniversalVariables(r0_vec,r_vec,tt,tm)
% This function solves the Lambert Problem using Universal Variables, as
% decribed by Vallado in Fundamentals of Astrodynamics and Applications.

% r0 - Initiatl position vector
% r  - Final position vector
% tt - Transfer time
% mu - Standard gravitational parameter
% tm - Transfer method
%      +1 is short way (< 180 deg, sometimes aka prograde)
%      -1 is long way (> 180 deg, sometimes aka retrograde)

%--------------------------------------------------------------------------

% Generate magnitudes of the position vectors
r0 = norm(r0_vec); r = norm(r_vec);

% Find cos(DeltaV) where DeltaV is the angle between the two positions.
cdv = dot(r0_vec,r_vec)/(r0*r);

% Find sin(DeltaV). tm determines the direction of travel.
sdv = tm*sqrt(1-cdv^2);

% Find A
A = tm*sqrt(r*r0*(1+cdv));
if A == 0.0
    fprintf('Orbit cannot be calculated\n')
end

% Initialize psi_n
psi_n = 0; %psi_0 = 0, set arbitrarily according to Vallado
c2 = 1/2;
c3 = 1/6;

psi_up = 4*pi^2; psi_low = -4*pi;
end