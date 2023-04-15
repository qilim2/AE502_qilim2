clear; close all; clc;


fprintf('Problem 2\n');
a = 1; % DU
e = 0.5;
i = 45; % degree
timespan = 100; % TU
mu = 1; % gravitational parameter is unity
n = sqrt(mu/a^3); % mean motion

i = deg2rad(i); % radians

L = n*a^2;
G = L*sqrt(1-e^2);
H = G*cos(i);

Ldot = 0; Gdot = 0; Hdot = 0;
ldot = 1/(2*L^2); gdot = 0; hdot = 0;

stepsize = 1e-3;
ts = 0:stepsize:timespan;

r_p = a*(1-e);
h = (r_p*mu*(1+e))^0.5; % Angular momentum
RAAN = 0;
AOP = 0;

% Note: the only value that changes is the mean anomaly => change in true
% anomaly

% Calculate eccentric anomaly from mean anomaly via Laguerre Method
M = 0;
E = Laguerre(M,e);

% Calculate true anomaly from eccentric anomaly
% Note: true anomaly can be calculated directly from mean anomaly, but
% approximation errors can occur when eccentricity is not small.
Beta = e/(1+sqrt(1-e^2));
TA = E + 2*atan2(Beta*sin(E),1-Beta*cos(E));
[r_vec,v_vec] = oe2rv(h,e,RAAN,deg2rad(i),AOP,TA,mu);
r_keep = zeros(length(ts),3);
v_keep = zeros(length(ts),3);
TA_keep = zeros(length(ts),1);
r_keep(1,:) = r_vec; v_keep(1,:) = v_vec;
for timestep = 1:length(ts)
    M = M + ldot*stepsize;
    E = Laguerre(M,e);
    TA = E + 2*atan2(Beta*sin(E),1-Beta*cos(E));
    TA_keep(timestep) = TA;
    [r_vec,v_vec] = oe2rv(h,e,RAAN,deg2rad(i),AOP,TA,mu);
    r_keep(timestep+1,:) = r_vec;
    v_keep(timestep+1,:) = r_vec;
end


plot3(r_keep(:,1),r_keep(:,2),r_keep(:,3))

