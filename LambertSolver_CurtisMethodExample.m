% ˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜
% Example_5_02
% ˜˜˜˜˜˜˜˜˜˜˜˜
%
% This program uses Algorithm 5.2 to solve Lambert’s problem
% for the data provided in Example 5.2.
%
% deg - factor for converting between degrees and radians
% pi - 3.1415926...
% mu - gravitational parameter (kmˆ3/sˆ2)
% r1, r2 - initial and final position vectors (km)
% dt - time between r1 and r2 (s)
% string - = 'pro' if the orbit is prograde
% = 'retro' if the orbit is retrograde
% v1, v2 - initial and final velocity vectors (km/s)
% coe - orbital elements [h e RA incl w TA a]
% where h = angular momentum (kmˆ2/s)
% e = eccentricity
% RA = right ascension of the ascending node
% (rad)
% incl = orbit inclination (rad)
% w = argument of perigee (rad)
% TA = true anomaly (rad)
% a = semimajor axis (km)
% TA1 - Initial true anomaly
% TA2 - Final true anomaly
% T - period of an elliptic orbit
%
% User M-functions required: lambert, coe_from_sv
% -----------------------------------------------------------
clear;
global mu
deg = pi/180;
mu = 1.327e11; %398600;
%...Input data from Example 5.2:
r1 = [ 5000 10000 2100]*1.496e+8;
r2 = [-14600 2500 7000]*1.496e+8;
dt = 86400*406;
grade = 0; % 0 for prograde, 1 for retrograde
%...

%...Algorithm 5.2:
%[v1, v2] = lambert(r1, r2, dt, 'pro');
%[v1, v2] = LambertSolver_CurtisMethod(r1, r2, dt, 0, mu);
[v1, v2, exitflag] = LambertSolver_ILBG(r1, r2, dt, 0, mu);

%...Algorithm 4.1 (using r1 and v1):
coe = rv2oe(r1, v1, mu);
%...Save the initial true anomaly:
TA1 = coe(6);
%...Algorithm 4.1 (using r2 and v2):
coe = rv2oe(r2, v2, mu);
%...Save the final true anomaly:
TA2 = coe(6);
%...Echo the input data and output the results to the command window:
fprintf('---------------------------------------------------')
fprintf('\n Example 5.2: Lambert''s Problem\n')
fprintf('\n\n Input data:\n');
fprintf('\n Gravitational parameter (kmˆ3/sˆ2) = %g\n', mu)
fprintf('\n r1 (km) = [%g %g %g]', ...
r1(1), r1(2), r1(3))
fprintf('\n r2 (km) = [%g %g %g]', ...
r2(1), r2(2), r2(3))
fprintf('\n Elapsed time (s) = %g', dt);
fprintf('\n\n Solution:\n')
fprintf('\n v1 (km/s) = [%g %g %g]', ...
v1(1), v1(2), v1(3))
fprintf('\n v2 (km/s) = [%g %g %g]', ...
v2(1), v2(2), v2(3))
fprintf('\n\n Orbital elements:')
fprintf('\n Angular momentum (kmˆ2/s) = %g', coe(1))
fprintf('\n Eccentricity = %g', coe(2))
fprintf('\n Inclination (deg) = %g', coe(4)/deg)
fprintf('\n RA of ascending node (deg) = %g', coe(3)/deg)
fprintf('\n Argument of perigee (deg) = %g', coe(5)/deg)
fprintf('\n True anomaly initial (deg) = %g', TA1/deg)
fprintf('\n True anomaly final (deg) = %g', TA2/deg)
fprintf('\n Semimajor axis (km) = %g', coe(7))
fprintf('\n Periapse radius (km) = %g', ...
coe(1)^2/mu/(1 + coe(2)))
if coe(2)<1
T = 2*pi/sqrt(mu)*coe(7)^1.5;
fprintf('\n Period:')
fprintf('\n Seconds = %g', T)
fprintf('\n Minutes = %g', T/60)
fprintf('\n Hours = %g', T/3600)
fprintf('\n Days = %g', T/24/3600)
end
fprintf('\n-----------------------------------------------\n')
% ˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜
