clear; close all; clc;
fprintf('AE 502 HW1 |\nQi Lim     |\n------------\n');
format long
fprintf('\nSOLVING LAMBERT''S PROBLEM\nThis code uses units of AU and days.\n\n');

% PROBLEM STATEMENT:
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
%-------------------------------------------------------------------------------

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
day2sec = 86400; %s/day
AU2km = 149597870.7; %km/AU
AU2m = AU2km * 1e3; %m/AU
muS = 1.32712440018e20; %m3/s2
%muS = muS / 1000^3; %km3/s2
muS = muS / AU2m^3 * day2sec^2; %Convert muS from m3/s2 to AU3/day2
%muS = muS * (6.68459e-12)^3 * (86400)^2;
deg = pi/180;

%% Convert to orbital elements (for problem 5)

r_vecs = [r1I, r2I, rE];
v_vecs = [v1I, v2I, vE];

for ii = 1:3
    
    coe = rv2oe(r_vecs(:,ii),v_vecs(:,ii),muS);
    
    if ii == 1
        fprintf('\n\n1I Orbital Elements:');
    elseif ii == 2
        fprintf('\n\n2I Orbital Elements:');
    elseif ii == 3
        fprintf('\n\nrE Orbital Elements:');
    else
        fprintf('INVALID INPUT. DISREGARD THE FOLLOWING\n\n');
    end

    fprintf('\n Angular momentum (AUˆ2/day) = %g', coe(1))
    fprintf('\n Eccentricity = %g', coe(2))
    fprintf('\n Inclination (deg) = %g', coe(4)/deg)
    fprintf('\n RA of ascending node (deg) = %g', coe(3)/deg)
    fprintf('\n Argument of perigee (deg) = %g', coe(5)/deg)
    fprintf('\n True anomaly initial (deg) = %g', coe(6)/deg)
    fprintf('\n Semimajor axis (AU) = %g', coe(7))
    fprintf('\n Periapse radius (AU) = %g', coe(1)^2/muS/(1 + coe(2)))
    if coe(2)<1
        T = 2*pi/sqrt(muS)*coe(7)^1.5*day2sec;
        fprintf('\n Period:')
        fprintf('\n Seconds = %g', T)
        fprintf('\n Minutes = %g', T/60)
        fprintf('\n Hours = %g', T/3600)
        fprintf('\n Days = %g', T/24/3600)
    end
    
end
%-------------------------------------------------------------------------------

%% Generate state vectors

% Julian Dates
jan_01_2017 = 2457754.500000;
dec_31_2017 = 2458118.500000;
aug_01_2017 = 2457966.500000;
jan_01_2019 = 2458484.500000;
jul_31_2020 = 2459061.500000;
jun_01_2019 = 2458635.500000;
jan_01_2022 = 2459580.500000;

% Set possible start and end times
dif1 = 1; %difference
dif2 = 1;
startdate = jan_01_2022;
t0_1Is = (jan_01_2017 - startdate):dif1:(dec_31_2017 - startdate); %day
tf_1Is = (aug_01_2017 - startdate):dif1:(jan_01_2019 - startdate); %day
t0_2Is = (jan_01_2017 - startdate):dif2:(jul_31_2020 - startdate); %day
tf_2Is = (jun_01_2019 - startdate):dif2:(jan_01_2022 - startdate); %day

% Produce state vectors of Earth for each departure date (km and km/s)
if length(t0_1Is) > length(t0_2Is)
    r0s = zeros(3,length(t0_1Is)); %Initialize matrix of possible r0s
    v0s = r0s; %Initialize matrix of possible v0s
    r0 = rE; v0 = vE;
    it1 = 1;
    while it1 <= length(t0_1Is)
        %fprintf('Iteration: %i\n',it1);
        [r_vec,v_vec] = UniversalOrbitPropagator(r0,v0,t0_1Is(it1),muS); %Propagate departure body to departure time
        r0s(:,it1) = r_vec;
        v0s(:,it1) = v_vec;
        r0 = r_vec; v0 = v_vec;
        it1 = it1 + 1;
    end
else
    r0s = zeros(3,length(t0_2Is)); %Initialize matrix of possible r0s
    v0s = r0s; %Initialize matrix of possible v0s
    r0 = rE; v0 = vE;
    it1 = 1;
    while it1 <= length(t0_2Is)
        %fprintf('Iteration: %i\n',it1);
        [r_vec,v_vec] = UniversalOrbitPropagator(r0,v0,t0_2Is(it1),muS); %Propagate departure body to departure time
        r0s(:,it1) = r_vec;
        v0s(:,it1) = v_vec;
        r0 = r_vec; v0 = v_vec;
        it1 = it1 + 1;
    end
end


fprintf('%g Departure State Vectors Calculated...\n', length(t0_1Is));

% Produce state vectors of 'Oumouamoua and Borisov for each arrival date (km and km/s)
rf1Is = zeros(3, length(tf_1Is)); %Initialize matrix of possible r11Is
vf1Is = rf1Is; %Initialize matrix of possible v11Is
rf2Is = zeros(3, length(tf_2Is)); %Initialize matrix of possible r12Is
vf2Is = rf2Is; %Initialize matrix of possible v12Is
rf1I = r1I; vf1I = v1I;
rf2I = r2I; vf2I = v2I;
it2 = 1;
while it2 <= length(tf_1Is) % 1I
  [r_vec,v_vec] = UniversalOrbitPropagator(rf1I,vf1I,tf_1Is(it2),muS); %Propagate 1I to arrival time
  rf1Is(:,it2) = r_vec;
  vf1Is(:,it2) = v_vec;
  rf1I = r_vec; vf1I = v_vec;
  it2 = it2 + 1;
end
fprintf('%g Arrival State Vectors Calculated for 1I...\n', length(tf_1Is));
while it2 <= length(tf_2Is) % 2I
  [r_vec,v_vec] = UniversalOrbitPropagator(rf2I,vf2I,tf_2Is(it2),muS); %Propagate 2I to arrival time
  rf2Is(:,it2) = r_vec;
  vf2Is(:,it2) = v_vec;
  rf2I = r_vec; vf2I = v_vec;
  it2 = it2 + 1;
end
fprintf('%g Arrival State Vectors Calculated for 2I...\n', length(tf_2Is));

%-------------------------------------------------------------------------------

%% Solving

warning('off','all')
%warning

% Solve for 1I
fprintf('\nSolving for 1I');
fprintf('\nt0 (out of %g) =  ', length(t0_1Is));

%Initialize matrices to hold total Delta-Vs
DepartureDeltaVs1I = zeros(length(t0_1Is),length(tf_1Is));
ArrivalDeltaVs1I = DepartureDeltaVs1I;
DeltaVsMatrix1I = DepartureDeltaVs1I; 
checks = DepartureDeltaVs1I;
for it3 = 1:length(t0_1Is)

    r1_DepartureBody = r0s(:,it3); %Radial position of departure body at departure
    v1_DepartureBody = v0s(:,it3); %Velocity of departure body (Earth in this case) at departure
    t0 = t0_1Is(it3); %Departure time

    % Print a nice counter
    if it3 <= 10
      fprintf('\b%i',it3);
    elseif it3 > 10 && it3 <= 100
      fprintf('\b\b%i',it3);
    elseif it3 > 100 && it3 <= 1000
      fprintf('\b\b\b%i',it3);
    else
      fprintf('\b\b\b\b%i',it3);
    end

    for it4 = 1:length(tf_1Is)

        tf = tf_1Is(it4); %Arrival time
        %fprintf('t0 = %g, tf = %g\n',t0,tf);
        TOF = abs(tf - t0); %Time of flight
        r2_ArrivalBody = rf1Is(:,it4); %Radial position of target body at arrival
        v2_ArrivalBody = vf1Is(:,it4); %Velocity of target body at arrival
        [v1_Spacecraft, v2_Spacecraft, check1] = LambertSolver_IzzoMethod(r1_DepartureBody, r2_ArrivalBody, TOF, 0, muS);
        %[v1_Spacecraft, v2_Spacecraft] = LambertSolver_CurtisMethod(r1_DepartureBody, r2_ArrivalBody, TOF, 0, muS);
            % v1_Spacecraft: velocity needed for transfer
            % v2_Spacecraft: velocity at arrival

        % Calculate total Delta-V
        DepartureDeltaV = abs(norm(v1_Spacecraft - v1_DepartureBody));
        ArrivalDeltaV = abs(norm(v2_ArrivalBody - v2_Spacecraft));
        TotalDeltaV = DepartureDeltaV + ArrivalDeltaV;

        % Convert total Delta-V from AU/day to km/s
        DepartureDeltaV = DepartureDeltaV * AU2km / day2sec;
        ArrivalDeltaV = ArrivalDeltaV * AU2km / day2sec;
        TotalDeltaV = TotalDeltaV * AU2km / day2sec;

        % Store in matrix
        DepartureDeltaVs1I(it3,it4) = DepartureDeltaV;
        ArrivalDeltaVs1I(it3,it4) = ArrivalDeltaV;
        DeltaVsMatrix1I(it3,it4) = TotalDeltaV;
        checks(it3,it4) = check1;
    end
end
fprintf('\nDelta-V Calculations for 1I Complete\n');

% Solve for 2I
fprintf('\nSolving for 2I');
fprintf('\nt0 (out of %g) =  ', length(t0_2Is));

%Initialize matrices to hold total Delta-Vs
DepartureDeltaVs2I = zeros(length(t0_2Is),length(tf_2Is));
ArrivalDeltaVs2I = DepartureDeltaVs2I;
DeltaVsMatrix2I = DepartureDeltaVs2I; 
for it5 = 1:length(t0_2Is)
    r1_DepartureBody = r0s(:,it5); %Radial position of departure body at departure
    v1_DepartureBody = v0s(:,it5); %Velocity of departure body (Earth in this case) at departure
    t0 = t0_2Is(it5); %Departure time

    % Print a nice counter
    if it5 <= 10
      fprintf('\b%i',it5);
    elseif it5 > 10 && it5 <= 100
      fprintf('\b\b%i',it5);
    elseif it5 > 100 && it5 <= 1000
      fprintf('\b\b\b%i',it5);
    else
      fprintf('\b\b\b\b%i',it5);
    end

    for it6 = 1:length(tf_2Is)

        tf = tf_2Is(it6); %Arrival time
        TOF = abs(tf - t0); %Time of flight
        r2_ArrivalBody = rf2Is(:,it6); %Radial position of target body at arrival
        v2_ArrivalBody = vf2Is(:,it6); %Velocity of target body at arrival
        [v1_Spacecraft, v2_Spacecraft] = LambertSolver_IzzoMethod(r1_DepartureBody, r2_ArrivalBody, TOF, 0, muS);
        % v1_Spacecraft: velocity needed for transfer
        % v2_Spacecraft: velocity at arrival

        % Calculate total Delta-V
        DepartureDeltaV = abs(norm(v1_Spacecraft - v1_DepartureBody));
        ArrivalDeltaV = abs(norm(v2_ArrivalBody - v2_Spacecraft));
        TotalDeltaV = DepartureDeltaV + ArrivalDeltaV;

        % Convert total Delta-V from AU/day to km/s
        DepartureDeltaV = DepartureDeltaV * AU2km / day2sec;
        ArrivalDeltaV = ArrivalDeltaV * AU2km / day2sec;
        TotalDeltaV = TotalDeltaV * AU2km / day2sec;

        % Store in matrix
        DepartureDeltaVs2I(it5,it6) = DepartureDeltaV;
        ArrivalDeltaVs2I(it5,it6) = ArrivalDeltaV;
        DeltaVsMatrix2I(it5,it6) = TotalDeltaV;
    end
end
fprintf('\nDelta-V Calculations for 2I Complete\n');

%-------------------------------------------------------------------------------

%% Finalization
fprintf('\n\nRUN COMPLETE\n\n');

% Plotting
plotting1 = DepartureDeltaVs1I;
plotting2 = DeltaVsMatrix1I;
[nr1I,nc1I] = size(DeltaVsMatrix1I);
for jj = 1:nr1I
    for kk = 1:nc1I
        orbitchecker = DeltaVsMatrix1I(jj,kk);
        if orbitchecker < 50 % Rendezvous
            plotting1(jj,kk) = orbitchecker;
        elseif orbitchecker < 20 % Flyby
            plotting2(jj,kk) = orbitchecker;
        else % Neither
            plotting1(jj,kk) = 0;
            plotting2(jj,kk) = 0;
        end
    end
end
%sum(plotting1,'all')
%sum(plotting2,'all')

figure(1)
contour(t0_1Is,tf_1Is,plotting1')
%xlim([t0_1Is(1),t0_1Is(end)]); ylim([tf_1Is(1),tf_1Is(end)]);
xlabel('Departure day'); ylabel('Arrival day'); zlabel('Delta-V');
figure(2)
contour(t0_1Is,tf_1Is,plotting2')
%xlim([t0_1Is(1),t0_1Is(end)]); ylim([tf_1Is(1),tf_1Is(end)]);
xlabel('Departure day'); ylabel('Arrival day'); zlabel('Delta-V');

plotting3 = DeltaVsMatrix2I;
plotting4 = DeltaVsMatrix2I;
[nr2I,nc2I] = size(DeltaVsMatrix2I);
for jj = 1:nr2I
    for kk = 1:nc2I
        orbitchecker = DeltaVsMatrix2I(jj,kk);
        if orbitchecker < 60 % Rendezvous
            plotting3(jj,kk) = orbitchecker;
        elseif orbitchecker < 20 % Flyby
            plotting4(jj,kk) = orbitchecker;
        else % Neither
            plotting3(jj,kk) = 0;
            plotting4(jj,kk) = 0;
        end
    end
end
%sum(plotting3,'all')
%sum(plotting4,'all')

figure(3)
contour(t0_2Is,tf_2Is,plotting3')
%xlim([t0_2Is(1),t0_2Is(end)]); ylim([tf_2Is(1),tf_2Is(end)]);
xlabel('Departure day'); ylabel('Arrival day'); zlabel('Delta-V');
figure(4)
contour(t0_2Is,tf_2Is,plotting4')
%xlim([t0_2Is(1),t0_2Is(end)]); ylim([tf_2Is(1),tf_2Is(end)]);
xlabel('Departure day'); ylabel('Arrival day'); zlabel('Delta-V');

