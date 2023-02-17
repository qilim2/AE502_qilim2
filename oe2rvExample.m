clear; clc;

deg = pi/180;
mu = 398600;
%...Input data (angles in degrees):
h = 80000;
e = 1.4;
RA = 40;
incl = 30;
w = 60;
TA = 30;
%...
%coe = [h, e, RA*deg, incl*deg, w*deg, TA*deg];
%...Algorithm 4.2 (requires angular elements be in radians):
[r, v] = oe2rv(h,e,RA*deg,incl*deg,w*deg,TA*deg,mu);
%...Echo the input data and output the results to the command window:
fprintf('---------------------------------------------------')
fprintf('\n Example 4.5\n')
fprintf('\n Gravitational parameter (kmˆ3/sˆ2) = %g\n', mu)
fprintf('\n Angular momentum (kmˆ2/s) = %g', h)
fprintf('\n Eccentricity = %g', e)
fprintf('\n Right ascension (deg) = %g', RA)
fprintf('\n Argument of perigee (deg) = %g', w)
fprintf('\n True anomaly (deg) = %g', TA)
fprintf('\n\n State vector:')
fprintf('\n r (km) = [%g %g %g]', r(1), r(2), r(3))
fprintf('\n v (km/s) = [%g %g %g]', v(1), v(2), v(3))
fprintf('\n-----------------------------------------------\n')
% ˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜
