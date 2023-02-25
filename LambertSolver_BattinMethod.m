function [outputArg1,outputArg2] = LambertSolver_BattinMethod(r0_vec,r_vec,tt,mu,tm)
% This function solves the Lambert Problem following Battin's method, as
% decribed by Vallado in Fundamentals of Astrodynamics and Applications.

% r0 - Initiatl position vector
% r  - Final position vector
% tt - Transfer time
% mu - Standard gravitational parameter
% tm - Transfer method
%      +1 is short way (< 180 deg, sometimes aka prograde)
%      -1 is long way (> 180 deg, sometimes aka retrograde)

%--------------------------------------------------------------------------

% INCOMPLETE


% Generate magnitudes of the position vectors
r0 = norm(r0_vec); r = norm(r_vec);

% Find cos(DeltaV) 
cdv = dot(r0_vec,r_vec)/(r0*r);

% Find DeltaV, is the angle between the two positions.
DeltaV = acos(cdv);

% Find sin(DeltaV). tm determines the direction of travel.
sdv = tm*sqrt(1-cdv^2);

% Find the chord length via the cosine law
c = sqrt(r0^2 + r2^2 - 2*r0*r*cdv);

% Define the semiperimeter, s (not the same as semiparameter)
s = (r0 + r + c)/2;

% Define epsilon
eps = (r - r0)/r0;

% Solve for tan^2(2w)
tan2w2 = (eps^2/4)/(sqrt(r/r0) + r/r0*(2 + sqrt(r/r0)));

% Find parabolic mean point radius
r_op = sqrt(r0*r)*(cos(DeltaV/4)^2 + tan2w2);


% Unlabeled ---------------------------------------------------------------
if tm == 1
    l = (sin(DeltaV/4)^2 + tan2w2)/(sin(DeltaV/4)^2 + tan2w2 + cos(DeltaV/2));
elseif tm == -1
    l = (cos(DeltaV/4)^2 + tan2w2 - cos(DeltaV/2))/(cos(DeltaV/4)^2 + tan2w2);
end

m = mu*tt^2/(8*r_op^3);
% Unlabeled ---------------------------------------------------------------


% Solve semimajor axis and semimparameter for minimum energy
a_min = s/2;
p_min = r0*r/c*(1-cdv);

% Solve for eccentricity of minimum-energy transfer
e_min = sqrt(1-2*p_min/s);

if e > 0
    x = l;
else
    x = 0.0;
end




% Generate alpha and beta
alpha_e = 2*asin(sqrt(s/(2*a_min)));
beta_e = 2*asin(sqrt((s-ec)/s));

% Find transfer time for minimizing semimajor axis
if tm == 1
    t_min_amin = sqrt(a_min^3/mu)*(alpha_e - (beta_e - sin(beta_e)));
elseif tm == -1
    t_min_amin = sqrt(a_min^3/mu)*(alpha_e + (beta_e - sin(beta_e)));
end

% Find absolute minimum transfer time
t_min_abs = (1/3)*sqrt(2/mu)*(s^1.5 - (s - c)^1.5);

% Eccentric Anomaly and Hyperbolic Anomaly
E = alpha_e - beta_e;
H = alpha_h - beta_h;





end