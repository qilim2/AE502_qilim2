
function J2_Perturbation_Encke
clear; close all; clc;
% This script follows Curtis' Example_12_02.m to solve a J2 perturbation of
% an orbit using Encke's method
% Conversions
hour2second = 60*60;
day2second = 24*hour2second;
degree2radian = pi/180; 
% I know that MATLAB has deg2rad as a built in function, but this should 
% keep it robust if I decide to port this to a different language.

% Constants
mu_Earth = 398600; %Gravitaional paramter in [km3/s2]
R_Earth = 6378; %Earth radius in [km]
J2_Earth = 1082.63e-6; %Earth J2 perturbation

% Initial orbital parameters (Example_12_02)
zp0 = 300; %Perigee altitude in [km]
za0 = 3062; %Apogee altitude in [km]
RA0 = 45*degree2radian; %Right ascension of the node in [rad]
i0 = 28*degree2radian; %Inclination in [rad]
w0 = 30*degree2radian; %Argument of perigee in [rad]
TA0 = 40*degree2radian; %True anomaly in [rad]

%{
% Given values for HW2
a = 26600; %semimajor-axis in km
i = 1.10654; %inclination in rad
e = 0.74; %eccentricity
omega = 5; %argument of perigee in deg
Omega = 90; %right ascension of the ascending node in deg
M0 = 10; %initial mean anomaly in deg
%}

% Inferred initial orbital paramters (Example_12_02)
rp0 = R_Earth + zp0; %Perigee radius in [km]
ra0 = R_Earth + za0; %Apogee radius in [km]
e0 = (ra0 - rp0)/(ra0 + rp0); %Eccentricity
a0 = (ra0 + rp0)/2; %Semimajor axis in [km]
h0 = sqrt(rp0*mu_Earth*(1+e0)); %Angular momentum in [km2/s]
T0 = 2*pi/sqrt(mu_Earth)*a0^1.5; %Period in [s]

t0 = 0; %Start time                                               --- INPUT
tf = 2*day2second; %End time                                      --- INPUT

% Store orbital elements and obtain initial state vector
coe0 = [h0 e0 RA0 i0 w0 TA0];
[r0_vec, v0_vec] = oe2rv(h0, e0, RA0, i0, w0, TA0, mu_Earth);
r0 = norm(r0_vec);
v0 = norm(v0_vec);

% Prepare for Encke procedure
del_t = T0/100; %Time step
options = odeset('MaxStep', del_t);

% Begin Encke integration
t = t0; %Initialize time scalar
tsave = t0; %Initialize vector of solution times
y = [r0_vec v0_vec]; %Initialize state vector (row vector)
del_y0 = zeros(6,1); %Initialize state vector perturbation

t = t + del_t; % First time step

% Loop over the time interval [t0, tf] with equal increments del_t
while t <= tf + del_t/2
    % Integrate eq. 12.7 over the time increment delt_t
    [dum,z] = ode45(@rates, [t0 t], del_y0, options);

    % Compute the osculating state vector at time t
    [rOsc_vec, vOsc_vec] = r0v02rv(r0_vec,v0_vec,t-t0,mu_Earth);

    % Rectify
    r0_vec = rOsc_vec + z(end,1:3);
    v0_vec = vOsc_vec + z(end,4:6);

    % Prepare for next time step
    tsave = [tsave; t];
    y = [y; [r0_vec v0_vec]];
    t = t + del_t;
    del_y0 = zeros(6,1);
end

t = tsave; %t is now the vector of equispaced solution times

% Extract the orbital elements from the state vector at each solution time
for it1 = 1:length(t)
    r_vec = [y(it1,1:3)];
    v_vec = [y(it1,4:6)];
    r(it1) = norm(r_vec);
    v(it1) = norm(v_vec);
    coe = rv2oe(r_vec,v_vec,mu_Earth);
    h(it1) = coe(1);
    e(it1) = coe(2);
    RA(it1) = coe(3);
    i(it1) = coe(4);
    w(it1) = coe(5);
    TA(it1) = coe(6);
end

%% Plotting
figure(1)
subplot(2,1,1)
plot(t/3600,(RA - RA0)/degree2radian)
title('Variation of Right Ascension')
xlabel('hours')
ylabel('{\it\Delta\Omega} (deg)')
grid on
grid minor
axis tight

subplot(2,1,2)
plot(t/3600,(w - w0)/degree2radian)
title('Variation of Argument of Perigee')
xlabel('hours')
ylabel('{\it\Delta\omega} (deg)')
grid on
grid minor
axis tight

figure(2)
subplot(3,1,1)
plot(t/3600,h - h0)
title('Variation of Angular Momentum')
xlabel('hours')
ylabel('{\it\Deltah} (km^2/s)')
grid on
grid minor
axis tight

subplot(3,1,2)
plot(t/3600,e - e0)
title('Variation of Eccentricity')
xlabel('hours')
ylabel('\it\Deltae')
grid on
grid minor
axis tight

subplot(3,1,3)
plot(t/3600,(i - i0)/degree2radian)
title('Variation of Inclination')
xlabel('hours')
ylabel('{\it\Deltai} (deg)')
grid on
grid minor
axis tight


%% Subfunctions
function dfdt = rates(t,f)
%wwwwwwwwwwwwwwwwwwwwwwwww
%
% This function calculates the time rates of Encke’s deviation in position
% del_r and velocity del_v.
% -----------------------------------------------------------------------

del_r = f(1:3)'; %Position deviation
del_v = f(4:6)'; %Velocity deviation

%...Compute the state vector on the osculating orbit at time t
% (Equation 12.5) using Algorithm 3.4:
[rOsc_vec,vOsc_vec] = r0v02rv(r0_vec, v0_vec, t-t0,mu_Earth);

%...Calculate the components of the state vector on the perturbed orbit
% and their magnitudes:
Rpp = rOsc_vec + del_r;
Vpp = vOsc_vec + del_v;
rosc = norm(rOsc_vec);
rpp = norm(Rpp);

%...Compute the J2 perturbing acceleration from Equation 12.30:
xx = Rpp(1); yy = Rpp(2); zz = Rpp(3);
fac = 3/2*J2_Earth*(mu_Earth/rpp^2)*(R_Earth/rpp)^2;
ap = -fac*[(1 - 5*(zz/rpp)^2)*(xx/rpp) ...
    (1 - 5*(zz/rpp)^2)*(yy/rpp) ...
    (3 - 5*(zz/rpp)^2)*(zz/rpp)];

%...Compute the total perturbing ecceleration from Equation 12.7:
F = 1 - (rosc/rpp)^3;
del_a = -mu_Earth/rosc^3*(del_r - F*Rpp) + ap;
dfdt = [del_v(1) del_v(2) del_v(3) del_a(1) del_a(2) del_a(3)]';
dfdt = [del_v del_a]'; %Return the derivative velocity and acceleration
%to ode45.

end %rates
end