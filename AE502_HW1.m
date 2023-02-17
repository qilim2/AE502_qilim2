clear; close all; clc;
fprintf('AE 502 HW1 |\nQi Lim     |\n------------\n');
format long
fprintf('\nSOLVING LAMBERT''S PROBLEM\nThis code uses units of AU and days.\n\n');
%{
Interstellar objects passing through our Solar System have been one of the most
amazing scientific discoveries of the last decade. Interstellar objects originate
beyond our Solar System and travel on hyperbolic trajectories with respect
to our Sun. 1I/’Oumouamoua (Figure 1) was the first interstellar object to
be discovered by Robert Weryk on October 19, 2017 [Meech et al., 2017]. In
2019, Comet 2I/Borisov was the second object discovered to be on a hyperbolic
trajectory through our Solar System [Opitom et al., 2019]. Getting a closer look
at those objects would yield invaluable scientific insights into how asteroids and
comets form around other suns. As future trajectory design engineers it now
falls to you to determine whether a rendez-vous or fly-by mission to such objects
is feasible! To this end you have to
1. (20 points) Write your own Universal Variable two-body orbit propagator
(e.g., Curtis’ Algorithms 3.3 and 3.4 [Curtis, 2013])!
2. (20 points) Write your own Lambert-Solver, based on e.g., Curtis’ Algorithm
5.2, or Gooding [1990] or even better Izzo [2015].
3. (20 points) Generate two pork chop plots [e.g. Bombardelli et al., 2018]
with departure dates from Earth in the range of January 2017 - December
2017 vs arrival dates in August 2017 - January 2019 vs total
ΔV color coded) for rendez-vous (ΔV < 50km/s) and fly-by missions
(ΔV < 20km/s) to 1I/’Oumouamoua, respectively.
4. (20 points) Generate two pork chop plots with departure dates from Earth
in the range of January 2017 - July 2020 vs arrival dates in June 2019 -
January 2022 vs total ΔV color coded) for rendez-vous (ΔV < 60km/s)
and fly-by missions (ΔV < 20km/s) to 2I/Borisov, respectively.
5. (10 points) Convert the initial state vectors to suitable orbital elements
and show that those objects are indeed interstellar!
6. (10 points) How realistic are those mission scenarios? If you had to pick
one, which one would it be?

Initial state vectors for the epoch of 2017-Jan-01 00:00:00.0000 UTC
i.e. JD 2457754.5., are in units of au and au/day:
r1I = [3.515868886595499 × 10−2,−3.162046390773074, 4.493983111703389]
v1I = [−2.317577766980901 × 10−3, 9.843360903693031 × 10−3,−1.541856855538041 × 10−2],
r2I = [7.249472033259724, 14.61063037906177, 14.24274452216359]
v2I = [−8.241709369476881 × 10−3,−1.156219024581502 × 10−2,−1.317135977481448 × 10−2],
and the initial state vector for the Earth is
rE = [−1.796136509111975 × 10−1, 9.667949206859814 × 10−1,−3.668681017942158 × 10−5]
vE = [−1.720038360888334 × 10−2,−3.211186197806460 × 10−3, 7.927736735960840 × 10−7].
Assume two-body motion with respect to the the Sun! Please upload your code
to your own GitHub repository! Good hunting!
%}

%% Initial State Vectors in AU and AU/day
% Initial state vectors for 'Oumouamoua
r1I = [3.515868886595499e-2;-3.162046390773074;4.493983111703389];
v1I = [-2.317577766980901e-3;9.843360903693031e-3;-1.541856855538041e-2];

% Initial state vectors for Borisov
r2I = [7.249472033259724;14.61063037906177;14.24274452216359];
v2I = [-8.241709369476881e-3;-1.156219024581502e-2;-1.317135977481448e-2];

% Initial state vectors for Earth
rE = [-1.796136509111975e-1; 9.667949206859814e-1;-3.668681017942158e-5];
vE = [-1.720038360888334e-2;-3.211186197806460e-3; 7.927736735960840e-7];

%% Constants
%day2sec = 86400; %s/day
%AU2km = 149597870.7; %km/AU
muS = 1.32712440018e20; %m3/s2
muS = muS * (6.68459e-12)^3 * (86400)^2;

%-------------------------------------------------------------------------------

%% Generate state vectors
% Set possible start times
dif = 1; %difference
t0 = 0:dif:365; %day

% Produce state vectors of Earth for each departure date (km and km/s)
r0s = zeros(3,length(t0)); %Initialize matrix of possible r0s
v0s = r0s; %Initialize matrix of possible v0s
r0 = rE; v0 = vE;
it1 = 1;

while it1 <= length(t0)
  %fprintf('Iteration: %i\n',it1);
  [r_vec,v_vec] = UniversalOrbitPropagator(r0,v0,t0(it1),muS);
  r0s(:,it1) = r_vec;
  v0s(:,it1) = v_vec;
  r0 = r_vec; v0 = v_vec;
  it1 = it1 + 1;
end
fprintf('Departure State Vectors Calculated...\n');

% Set possible end times
tf = 213:dif:365*2;

% Produce state vectors of 'Oumouamoua and Borisov for each arrival date (km and km/s)
rf1Is = zeros(3, length(tf)); %Initialize matrix of possible r11Is
vf1Is = rf1Is; %Initialize matrix of possible v11Is
rf2Is = zeros(3, length(tf)); %Initialize matrix of possible r12Is
vf2Is = rf2Is; %Initialize matrix of possible v12Is
rf1I = r1I; vf1I = v1I;
rf2I = r2I; vf2I = v2I;
it2 = 1;

while it2 <= length(tf)
  [r_vec,v_vec] = UniversalOrbitPropagator(rf1I,vf1I,tf(it2),muS);
  rf1Is(:,it2) = r_vec;
  vf1Is(:,it2) = v_vec;
  rf1I = r_vec; vf1I = v_vec;
  [r_vec,v_vec] = UniversalOrbitPropagator(rf2I,vf2I,tf(it2),muS);
  rf2Is(:,it2) = r_vec;
  vf2Is(:,it2) = v_vec;
  rf2I = r_vec; vf2I = v_vec;
  it2 = it2 + 1;
end
fprintf('Arrival State Vectors Calculated...\n');

%-------------------------------------------------------------------------------

%% Solving





%-------------------------------------------------------------------------------

%% Finalization
fprintf('\n\nRUN COMPLETE\n\n');
