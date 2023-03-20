clear; close all; clc;
fprintf('AE 502 HW2 |\nQi Lim     |\n------------\n');
format long

% Problem Statement
%{
So-called “Molniya orbits” were discovered by Soviet scientists in the 1960s as an
alternative to geostationary orbits, which, when launched from high latitudes,
require large launch energies to achieve a high perigee and to change inclination
in order to orbit over the equator. A long duration stay over a target communications 
area could be achieved using highly elliptical orbits with an apogee over
the desired territory, see Figure 1. The orbit’s name refers to the “lightning”
speed with which the satellite passes through perigee. In order to minimize
station keeping fuel expenditure Molniya orbits were designed to be “frozen”
under Earth s J2 perturbation and keep the nodal precession small.
%}

% Variables
%{
a = semimajor axis [km]
e = eccentricity
incl = inclination
n_prec = nodal precession rate
J2_Earth = J2 perturbation for Earth
R_Earth = Earth radius [km]
GM_earth = standard gravitational parameter of Earth [km3/s2]
%}


%% Constants
J2_Earth = 0.00108;
R_Earth = 6370; %km
GM_Earth = 3.986e5; %km3/s2
J2_Mars = 0.00196;
R_Mars = 3390; %km
GM_Mars = 4.282e4; %km3/s2

%% Problem 1
fprintf('\nProblem 1\n');
% Goal: find orbital elements a, e, and incl for a Molniya orbit and the
% lowest nodal precession rate.

% We want the satellite to orbit the Earth three times per day:
% 24 h * 60 min * 60 s / 3 orbits
Period_1 = 24*60*60/3; %s

% Minimum perigee altitude (given)
alt_p1_min = 600; %km
r_p1_min = alt_p1_min + R_Earth; %km

% Finding semimajor-axis from period equation
a_1 = (GM_Earth*(Period_1/(2*pi))^2)^(1/3); %km

% Finding eccentricity using semimajor-axis
e_1 = 1-r_p1_min/a_1;

% Calculate mean motion
n_1 = 2*pi/Period_1; %rad/s

% Find inclination from equation (1)
% Since n, J2, R, a, and e are non-zero, use only the numerator and solve
% for incl.
incl_p = acos(1/sqrt(5)); %rad
incl_n = acos(-1/sqrt(5)); %rad

% Nodal precession
n_prec_11 = -1.5*n_1*J2_Earth*(R_Earth/a_1)^2*cos(incl_p)/((1-e_1^2)^2); %rad/s
n_prec_12 = -1.5*n_1*J2_Earth*(R_Earth/a_1)^2*cos(incl_n)/((1-e_1^2)^2); %rad/s

fprintf('Semimajor-axis, a = %f km\n',a_1);
fprintf('Eccentricity, e = %f\n',e_1);
fprintf('Inclination 1, incl = %f rad = %f deg\n',incl_p,rad2deg(incl_p));
fprintf('Inclination 2, incl = %f rad = %f deg\n',incl_n,rad2deg(incl_n));
fprintf('Nodal Precession 1, n_prec = %e rad/s = %e deg/s\n',n_prec_11,rad2deg(n_prec_11));
fprintf('Nodal Precession 2, n_prec = %e rad/s = %e deg/s\n',n_prec_12,rad2deg(n_prec_12));

%% Problem 2
fprintf('\nProblem 2\n');
% Goal: find orbital elements, a, e, and incl for a Molniya orbit around
% Mars as well as the lowest nodal precession rate

% Orbit Mars once per Martian day
Period_2 = 24*60*60 + 39*60 + 35; %s

% Minimum periapse altitude (given)
alt_p2_min = 400; %km
r_p2_min = alt_p2_min + R_Mars; %km

% Finding semimajor-axis from period equation
a_2 = (GM_Mars*(Period_2/(2*pi))^2)^(1/3); %km

% Finding eccentricity using semimajor-axis
e_2 = 1-r_p2_min/a_2;

% Calculate mean motion
n_2 = 2*pi/Period_2; %rad/s

% Find inclination from equation (1)
% Since n, J2, R, a, and e are non-zero, use only the numerator and solve
% for incl.
% Same as for problem 1. Retain same values.

% Nodal precession
n_prec_21 = -1.5*n_2*J2_Mars*(R_Mars/a_2)^2*cos(incl_p)/((1-e_2^2)^2); %rad/s
n_prec_22 = -1.5*n_2*J2_Mars*(R_Mars/a_2)^2*cos(incl_n)/((1-e_2^2)^2); %rad/s

fprintf('Semimajor-axis, a = %f km\n',a_2);
fprintf('Eccentricity, e = %f\n',e_2);
fprintf('Inclination 1, incl = %f rad = %f deg\n',incl_p,rad2deg(incl_p));
fprintf('Inclination 2, incl = %f rad = %f deg\n',incl_n,rad2deg(incl_n));
fprintf('Nodal Precession 1, n_prec = %e rad/s = %e deg/s\n',n_prec_21,rad2deg(n_prec_21));
fprintf('Nodal Precession 2, n_prec = %e rad/s = %e deg/s\n',n_prec_22,rad2deg(n_prec_22));

%% Problem 3
fprintf('\nProblem 3\n');
% Goal: code up non-averaged perturbed equations of motion, either
% Planetary Equations or special perturbation equations, and solve
% numerically to study the effect of Earth's J2 perturbation on a Molniya
% orbit. Plot the time evolution of all orbital elements except the mean
% anomaly (M) for 100 days, then describe the results.

% Given values
a = 26600; %semimajor-axis in km
i = 1.10654; %inclination in rad
e = 0.74; %eccentricity
omega = 5; %argument of perigee in deg
Omega = 90; %right ascension of the ascending node in deg
M0 = 10; %initial mean anomaly in deg

% Convert values for input into function: J2_Perturbation_Gauss
r_p0 = a*(1-e); % Periapse radius
r_a0 = a*(1+e); % Apoapse radius
i0 = rad2deg(i);
omega0 = omega;
Omega0 = Omega;
R_body = R_Earth;
mu = GM_Earth;
J2 = J2_Earth;
tf = 100;

% Calculate eccentric anomaly from mean anomaly via Laguerre Method
E = Laguerre(deg2rad(M0),e);

% Calculate true anomaly from eccentric anomaly
% Note: true anomaly can be calculated directly from mean anomaly, but
% approximation errors can occur when eccentricity is not small.
Beta = e/(1+sqrt(1-e^2));
TA = E + 2*atan2(Beta*sin(E),1-Beta*cos(E));
TA0 = TA;

% Solve
J2_Perturbation_Gauss1(r_p0,r_a0,i0,omega0,Omega0,TA0,mu,J2,R_body,tf)
%J2_Perturbation_Gauss2(a,e,i0,omega0,Omega0,TA0,mu,J2,R_body,tf)

%{
% Test values (from Curtis Example)
mu = 398600;
RE = 6378;
J2 = 1082.63e-6;
rp0 = RE + 300;
ra0 = RE + 3062;
RA0 = 45;
i0 = 28;
w0 = 30;
TA0 = 40;
tf = 2;

J2_Perturbation_Gauss1(rp0,ra0,i0,w0,RA0,TA0,mu,J2,RE,tf)
%}

function E = Laguerre(M,e)
% This function solves the Eccentric Anomaly given the Mean Anomaly and the
% eccentricity of the orbit using Laguerre's method.

% See AIAA-86-0084, An Improved Algorithm Due to Laguerre for the Solution
% of Kepler's Equation, Conway, 1986
% and also, Prussing and Conway Orbital Mechanics, p. 38

% Starting value
u = M + e;
n  = 5;
E  = (M*(1-sin(u)) + u*sin(M)) / (1 + sin(M) - sin(u));

% Iterate
for it1 = 1:5
    num = n*(E-e*sin(E));
    den = (1-e*cos(E)) + sign(1-e*cos(E))*(abs((n-1)^2*(1-e*cos(E))^2 ...
        - n*(n-1)*(E-e*sin(E))*(e*sin(E))))^0.5;
    E = num/den;
end
end

