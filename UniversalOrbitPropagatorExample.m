% ˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜
% Example_3_07
% ˜˜˜˜˜˜˜˜˜˜˜˜
%
% This program computes the state vector (R,V) from the
% initial state vector (R0,V0) and the elapsed time using the
% data in Example 3.7.
%
% mu - gravitational parameter (e.g. kmˆ3/sˆ2)
% R0 - the initial position vector (e.g. km)
% V0 - the initial velocity vector (e.g. km/s)
% R - the final position vector (e.g. km)
% V - the final velocity vector (e.g. km/s)
% t - elapsed time (e.g. s)
%
% User M-functions required: rv_from_r0v0
% ------------------------------------------------------------
clear
%global mu
%mu = 398600;
%...Input data for Example 3.7:
%R0 = [10000 0 0];
%V0 = [0 9.2 0];
t = 1;


R0 = [-1.796136509111975e-1; 9.667949206859814e-1;-3.668681017942158e-5];
V0 = [-1.720038360888334e-2;-3.211186197806460e-3; 7.927736735960840e-7];
% Constants
muS = 1.32712440018e20; %m3/s2
muS = muS * (6.68459e-12)^3 * (86400)^2;

[R V] = UniversalOrbitPropagator(R0, V0, t, muS);

%...
%...Algorithm 3.4:
%[R V] = UniversalOrbitPropagator(R0, V0, t, mu);
%...Echo the input data and output the results to the command window:
fprintf('---------------------------------------------------')
fprintf('\n Example 3.7\n')
fprintf('\nmu = %f\n',muS);
fprintf('\n Initial position vector (km):')
fprintf('\n r0 = (%g, %g, %g)\n', R0(1), R0(2), R0(3))
fprintf('\n Initial velocity vector (km/s):')
fprintf('\n v0 = (%g, %g, %g)', V0(1), V0(2), V0(3))
fprintf('\n\n Elapsed time = %g s\n',t)
fprintf('\n Final position vector (km):')
fprintf('\n r = (%g, %g, %g)\n', R(1), R(2), R(3))
fprintf('\n Final velocity vector (km/s):')
fprintf('\n v = (%g, %g, %g)', V(1), V(2), V(3))
fprintf('\n-----------------------------------------------\n')
% ˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜

