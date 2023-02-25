clear; clc;

%% Constants
day2sec = 86400; %s/day
AU2km = 149597870.7; %km/AU
AU2m = AU2km * 1e3; %m/AU
%muS = 1.32712440018e20; %m3/s2
%muS = muS / AU2m^3 * day2sec^2; %Convert muS from m3/s2 to AU3/day2


%% Test values #1
fprintf('\nTest case 1\n');
r1 = [5644, 2830, 4170]; %km
r2 = [-2240, 7320, 4980]; %km
mu = 3.986004418e5; %km^3/s^2;
tof=1200; %s

fprintf('r1 = %s km\n',mat2str(r1));
fprintf('r2 = %s km\n',mat2str(r2));
fprintf('mu = %g km^3/s^2\ntof = %g seconds\n',mu,tof);

[v1, v2] = LambertSolver_IzzoMethod(r1,r2,tof,0,mu);
fprintf('v1 = %s km/s\n',mat2str(v1));
fprintf('v2 = %s km/s\n',mat2str(v2));

[v1, v2] = LambertSolver_CurtisMethod(r1,r2,tof,0,mu);
fprintf('v1 = %s km/s\n',mat2str(v1));
fprintf('v2 = %s km/s\n',mat2str(v2));



%% Test values #2
fprintf('\nTest case 2\n');
rE = [-1.796136509111975e-1; 9.667949206859814e-1;-3.668681017942158e-5];
r1I = [3.515868886595499e-2;-3.162046390773074;4.493983111703389];
v1I = [-2.317577766980901e-3;9.843360903693031e-3;-1.541856855538041e-2];
day2sec = 86400; %s/day
AU2km = 149597870.7; %km/AU
AU2m = AU2km * 1e3; %m/AU
muS = 1.32712440018e20; %m3/s2
muS = muS / AU2m^3 * day2sec^2; %Convert muS from m3/s2 to AU3/day2
tof = 400;
[r2,v2] = UniversalOrbitPropagator(r1I,v1I,tof,muS);
r1 = rE;

fprintf('r1 = %s AU\n',mat2str(r1));
fprintf('r2 = %s AU\n',mat2str(r2));
fprintf('mu = %g AU^3/day^2\ntof = %g days\n',muS,tof);

[v1, v2] = LambertSolver_IzzoMethod(r1,r2,tof,0,muS);
v1 = v1 * AU2km / day2sec;
v2 = v2 * AU2km / day2sec;
fprintf('v1 = %s km/s\n',mat2str(v1));
fprintf('v2 = %s km/s\n',mat2str(v2));

[v1, v2] = LambertSolver_CurtisMethod(r1,r2,tof,1,muS);
v1 = v1 * AU2km / day2sec;
v2 = v2 * AU2km / day2sec;
fprintf('v1 = %s km/s\n',mat2str(v1));
fprintf('v2 = %s km/s\n',mat2str(v2));