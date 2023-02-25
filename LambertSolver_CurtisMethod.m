function [v1_vec,v2_vec] = LambertSolver_CurtisMethod(r1_vec,r2_vec,t,grade,mu)

% grade: 0 for prograde, 1 for retrograde
% Requires function stumpff.m

% Inputs:
% r1_vec, initial position vector (km)
% r2_vec, final position vector (km)
% t, transfer time / time of flight (s)
% mu, gravitational parameter (km^3/s^2)

% Outputs:
% v1_vec, initial velocity vector (not including velocity of departure body,
% e.g. Earth orbit velocity around the Sun).
% v2_vec, final velocity vector (not including velocity of arrival body,
% e.g. asteroid velocity)


% This code solves Lambert's problem using the method described by Curtis.
% Transfer energy is minimized for the time of flight (t) given as an input.
% Total Delta-V is calculated as:
% norm(abs(v1)vec - vDepartureBody_vec)) + norm(abs(v2_vec - vArrivalBody_vec))

%------------------------------------------------------------------------------

% Ensure that r1_vec and r2_vec are column vectors
[r1_rows,r1_cols] = size(r1_vec);
[r2_rows,r2_cols] = size(r2_vec);
if r1_rows < r1_cols
    r1_vec = transpose(r1_vec);
end
if r2_rows < r2_cols
    r2_vec = transpose(r2_vec);
end

% Calculate norms
r1 = norm(r1_vec); r2 = norm(r2_vec);

% Change in true anomaly, theta
costheta = dot(r1_vec,r2_vec)/(r1*r2);
% If costheta > 0, then deltatheta lies in first or fourth quadrant.
% If costheta < 0, then deltatheta lies in second or third quadrant

%% Resolving Ambiguity
c12 = cross(r1_vec,r2_vec);
theta = acos(costheta);

% Modify theta based on prograde or retrograde trajectory, as well as the inclination.
if grade == 0
    if c12(3) < 0
        theta = 2*pi - theta;
    end
elseif grade == 1
    if c12(3) >= 0
        theta = 2*pi - theta;
    end
else
    fprintf('\n Prograde trajectory assumed.\n');
end

%% Solving

A = sin(theta)*sqrt(r1*r2/(1-costheta));
%disp(A)

% Find z through iteration
z = -100;
[S,C] = stumpff(z);
y = r1 + r2 + A*(z*S - 1)/sqrt(C);
%disp(y)
F = (y/C)^1.5*S + A*sqrt(y) - sqrt(mu)*t;
%disp(F)

while real(F) < 0
    z = z + 0.1;
    [S,C] = stumpff(z);
    y = r1 + r2 + A*(z*S - 1)/sqrt(C);
    F = (y/C)^1.5*S + A*sqrt(y) - sqrt(mu)*t;
    %fprintf('z iterate: %f\n',z);
end
%disp(F)

tol = 1.e-8; %tolerance
nmax = 10000;
ratio = 1;
n = 0;
while abs(ratio) > tol && n <= nmax
    n = n + 1;
    [S,C] = stumpff(z);
    y = (r1 + r2 + A*(z*S - 1)/sqrt(C));
    F = (y/C)^1.5*S + A*sqrt(y) - sqrt(mu)*t;
    if z == 0
        dFdz = sqrt(2)/40*y^1.5 + A/8*(sqrt(y) + A*sqrt(1/2/y));
    else
        dFdz = (y/C)^1.5*(1/2/z*(C - 3*S/2/C) + 3*S^2/4/C) ...
                + A/8*(3*S/C*sqrt(y) + A*sqrt(C/y));
    end
    ratio = F/dFdz;
    z = z - ratio;
end
%disp(z)
if n >= nmax
    fprintf('Max iterations exceeded. %g \n', nmax);
end

[S,C] = stumpff(z);
y = (r1 + r2 + A*(z*S - 1)/sqrt(C));
f = 1 - y/r1;
g  = A*sqrt(y/mu);
gdot = 1- y/r2;
v1_vec = 1/g*(r2_vec - f*r1_vec);
v2_vec = 1/g*(gdot*r2_vec - r1_vec);


end
