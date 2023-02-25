clear; clc;
% Example for UniversalKepler.m

mu = 398600; %km^3/s^2
r0 = 10000;
vr0 = 3.0752;
DeltaT = 3600;
a = -19655; %semimajor axis in km

chi = UniversalKepler(mu,DeltaT,r0,vr0,1/a);
