function J2_Perturbation_Gauss1(r_p0,r_a0,i0,omega0,Omega0,TA0,mu,J2,R_body,tf)
% This function solves J2 perturbations for an orbit over a given time
% period using MATLAB's ode45. This code is based on Curtis' Example_12_06
% with modifications.

%clear; close all; clc;
fprintf('Running function: J2_Perturbation_Gauss\n\n')

fprintf('Inputs:   (ensure units match)\n');
fprintf('r_p0 (periapse radius in [km] at t=0) = %f \n',r_p0);
fprintf('r_a0 (apoapse radius in [km] at t=0) = %f \n',r_a0);
fprintf('i0 (inclination in [deg] at t=0) = %f \n',i0);
fprintf('omega0 (Argument of Periapse in [deg] at t=0) = %f \n', omega0);
fprintf('Omega0 (Right Ascension of the Ascending Node in [deg] at t=0) = %f \n', Omega0);
fprintf('TA0 (True Anomaly in [deg] at t=0) = %f \n', TA0);
fprintf('mu (gravitational parameter (GM) in [km3/s2]) = %f \n',mu);
fprintf('J2 (J2 perturbation factor) = %f \n', J2);
fprintf('R_body (radius of the central body in [km]) = %f \n',R_body);
fprintf('tf (elapsed/final time in [days]) = %f \n', tf);

%% Preparation
% Conversion factors
% Utilize conversion factors to minimize use of MATLABS built-in functions.
hour2second = 60*60; % Multiply by this to convert from hours to seconds.
day2second = 24*hour2second; % Multiply by this to convert from days to seconds.
degree2radian = pi/180; % Multiply by this to convert from degrees to radians.

% Convert values
i0 = i0*degree2radian;
omega0 = omega0*degree2radian;
Omega0 = Omega0*degree2radian;
TA0 = TA0*degree2radian;

% Inferred orbital parameters
e0 = (r_a0 - r_p0)/(r_a0 + r_p0); % Eccentricity
h0 = (r_p0*mu*(1+e0))^0.5; % Angular momentum in [km2/s]
a0 = (r_p0 + r_a0)/2; % Semimajor axis in [km]
Period0 = 2*pi/(mu^0.5)*a0^1.5; % Period in [s]

% Store initial orbital elements
coe0 = [h0 e0 Omega0 i0 omega0 TA0];

%% Solve
% Integrate the Gauss variational equations from t0 to tf
t0 = 0;
tf = tf*day2second; % Convert final time to seconds
soln_points = 5000; % Number of solution points
tspan = linspace(t0,tf,soln_points);
options = odeset('reltol',1.e-8,'abstol',1.e-8,'initialstep',Period0/1000);
y0 = coe0';
[t,y] = ode45(@rates,tspan,y0,options);

% Save values for each variable
h = y(:,1);
e = y(:,2);
Omega = y(:,3);
i = y(:,4);
omega = y(:,5);
TA = y(:,6);

%% Plotting
figure(2)
subplot(5,1,1)
plot(t/3600,(Omega)/degree2radian)
title('Right Ascension of the Ascending Node (degrees)')
xlabel('hours')
grid on
grid minor
axis tight

subplot(5,1,2)
plot(t/3600,(omega)/degree2radian)
title('Argument of Periapse (degrees)')
xlabel('hours')
grid on
grid minor
axis tight

subplot(5,1,3)
plot(t/3600,h)
title('Angular Momentum (km^2/s)')
xlabel('hours')
grid on
grid minor
axis tight

subplot(5,1,4)
plot(t/3600,e)
title('Eccentricity')
xlabel('hours')
grid on
grid minor
axis tight

subplot(5,1,5)
plot(t/3600,(i)/degree2radian)
title('Inclination (degrees)')
xlabel('hours')
grid on
grid minor
axis tight

%% Subfunctions
    function dfdt = rates(t,f)
    % This function calculates the time rates of the orbital elements from
    % Gauss's variational equations. This code is primarily from Curtis'
    % example code. This code is checked against Curtis' explanations in
    % section 12.6 of Orbital Mechanics for Engineers, 2014.

    % Orbital elements at time t
    h = f(1); % Angular momentum
    e = f(2); % Eccentricity
    Omega = f(3); % Right ascension of the ascending node
    i = f(4); % Inclination
    omega = f(5); % Argument of periapse
    TA = f(6); % True anomaly

    % Inferred values
    r = h^2/mu/(1+e*cos(TA)); % Radius
    u = omega + TA; % Argument of latitude

    % Orbital element rates at time t
    hdot = -3/2*J2*mu*R_body^2/r^3*sin(i)^2*sin(2*u);

    edot = 3/2*J2*mu*R_body^2/h/r^3*(h^2/mu/r*sin(TA)*(3*sin(i)^2*sin(u)^2 - 1) ...
        - sin(2*u)*sin(i)^2*((2+e*cos(TA))*cos(TA)+e));

    TAdot =  h/r^2 + 3/2*J2*mu*R_body^2/e/h/r^3*(h^2/mu/r*cos(TA)*(3*sin(i)^2*sin(u)^2 - 1) ...
        + sin(2*u)*sin(i)^2*sin(TA)*(h^2/mu/r + 1));

    Omegadot = -3*J2*mu*R_body^2/h/r^3*sin(u)^2*cos(i);

    idot = -3/4*J2*mu*R_body^2/h/r^3*sin(2*u)*sin(2*i);

    omegadot = 3/2*J2*mu*R_body^2/e/h/r^3*(-h^2/mu/r*cos(TA)*(3*sin(i)^2*sin(u)^2 - 1) ...
        - sin(2*u)*sin(i)^2*sin(TA)*(2 + e*cos(TA)) + 2*e*cos(i)^2*sin(u)^2);

    % Pass back to ode45
    dfdt = [hdot edot Omegadot idot omegadot TAdot]';

    end
end