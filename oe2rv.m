function [r_vec,v_vec] = oe2rv(h,e,Omega,i,omega,TA,mu)
% This function is straight from Curtis
% All angles (Omega, i, omega) must be given in radians.

rp = (h^2/mu) * (1/(1 + e*cos(TA))) * (cos(TA)*[1;0;0] + sin(TA)*[0;1;0]);
vp = (mu/h) * (-sin(TA)*[1;0;0] + (e + cos(TA))*[0;1;0]);

%...Equation 4.39:
R3_W = [ cos(Omega) sin(Omega) 0
-sin(Omega) cos(Omega) 0
0 0 1];
%...Equation 4.40:
R1_i = [1 0 0
0 cos(i) sin(i)
0 -sin(i) cos(i)];
%...Equation 4.41:
R3_w = [ cos(omega) sin(omega) 0
-sin(omega) cos(omega) 0
0 0 1];
%...Equation 4.44:
Q_pX = R3_W'*R1_i'*R3_w';
%...Equations 4.46 (r and v are column vectors):
r_vec = (Q_pX*rp)';
v_vec = (Q_pX*vp)';


end