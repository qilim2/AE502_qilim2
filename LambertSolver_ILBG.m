function [v1_vec, v2_vec] = LambertSolver_ILBG(r1_vec, r2_vec, TOF, nOrbits, mu)

%% Statement
% This solver is based off of Rody Oldenhuis' code, which utilizes Izzo's
% method and Lancaster & Blanchard's method, with modifications by Gooding.

% I have rewritten this in order to better understand the full process.
% Variables have be renamed per my own understanding of the content. If a
% variable name remains consistent with the original code, this is because
% the name makes the most sense. I also make use of Vallado's Fundamentals
% of Astrodynamics and Applications 4th ed., and Curtis' Orbital Mechanics
% for Engineering Students for cross-reference.

% Name: Qi Lim
% Date: February, 2023
% Distribution: For use in AE 502, Advanced Orbital Mechanics .

% University of Illinois Urbana Champaign.

%% Citation of Oldenhuis' code
% Rody Oldenhuis, orcid.org/0000-0002-3162-3660. "Lambert" <version>,
% <date you last used it>. MATLAB Robust solver for Lambert's
% orbital-boundary value problem.
% https://nl.mathworks.com/matlabcentral/fileexchange/26348

%% Inputs and Outputs
% r1_vec: Initial position vector in [km] or [AU].
%   Size: [1x3] or [3x1]
% r2_vec: Final position vector in [km] or [AU].
%   Size: [1x3] or [3x1]
% nOrbits: Number of orbits to complete.
%   Size: [1x1]
% TOF: Time of flight in [days].
%   Size: [1x1]
% mu: Standard Gravitational Parameter of the central body in [km3/s2].
%   Size: [1x1]

% v1_vec: Initial velocity vector in [km/s].
%   Size: [1x3] or [3x1]
% v2_vec: Final velocity vector in [km/s].
%   Size: [1x3] or [3x1]

%% Preparation of position vectors
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

% Calculate unit vectors
%r1_vecUnit = r1_vec/r1;
r2_vecUnit = r2_vec/r2;

%% Initial values
tolerance = 1e-14;
check1 = false;
day2sec = 86400; %seconds/day

% Nondimensionalize by normalizing to initial vector.
r1_vecND = r1_vec/r1; %[ND]
r2_vecND = r2_vec/r1; %[ND]
r2_ND = norm(r2_vecND); %[ND]
V = sqrt(mu/r1); %orbital velocity [km/s]
T = r1/V; %time [s]
TOF = TOF*day2sec/T; %translate Time of Flight to seconds, then nondimensionalize.
%   May no longer be positive.

% Ensure within bounds
dtheta = acos( max(-1, min(1, dot(r1_vec,r2_vec)/r2_ND)));
% Note: Appears to be finding the transfer angle. However, this differs
% from what I find in Vallado and Curtis' books. That would be:
% dtheta = dot(r1vec,r2_vec)/(r1*r2)

% Branching
left_branch = sign(nOrbits); % Get sign
long_way = sign(TOF); % Get sign
nOrbits = abs(nOrbits);
TOF = abs(TOF);
if long_way < 0
    dtheta = 2*pi - dtheta;
end

% Note: it appears that most of the following shares similarities with
% Battin's method, as seen in Vallado.

%% Derived quantities
% Find chord length using the cosine law between the two position vectors
c = sqrt(1 + r2_ND^2 - 2*r2_ND*cos(dtheta)); %[ND]
    % This is just a version of c = sqrt(r1^2 + r2^2 - 2*r1*r2*cos(dtheta)^2)
    % but modified for non dimensionalized values.

% Define the semiperimeter
s = (1 + r2_ND + c)/2; %[ND]
    % Again, this is just a modified version from what is in Vallado's
    % book. He in turn references Prussing and Conway, as well as Kaplan.

% Semimajor axis for Minimum Energy transfer
a_min = s/2; %[ND]

% Lambda paramter (apparently from Battin)
lambda = sqrt(r2_ND)*cos(dtheta/2)/s;
    % This can be found on page 5 of Izzo's paper. Again, modified for ND.

% Nondimensionalize TOF
TOF_ND = sqrt(2*mu/(s^3))*TOF;

% Cross-product
cr1r2_vec = cross(r1_vecND, r2_vecND);
cr1r2 = norm(cr1r2_vec); %magnitude
cr1r2_uv = cr1r2_vec/cr1r2; %unit vector

%% Initial estimates

logTOF = log(TOF); % Compute logarithm of time of flight

if nOrbits == 0 % Single revolution
    % Initial values. These come from page 11 of Izzo's paper. The
    % hardcoded numbers come from the solution using non dimensionalized
    % values.
    iv1 = -0.5233; % First initial guess
    iv2 = +0.5233; % Second initial guess
    x1 = log(1+iv1); % Transformed first initial guess
    x2 = log(1+iv2); % Transformed second initial guess

else % Multiple revolutiosn (0, 1, or 2 solutions)
    if left_branch < 0
        iv1 = -0.5234; % First initial guess on left branch
        iv2 = +0.2234; % Second initial guess on left branch
    else
        iv1 = +0.7234; % First initial guess on right branch
        iv2 = +0.5234; % Second initial guess on right branch
    end

    x1 = tan(iv1*pi/2); % Transformed first initial guess
    x2 = tan(iv2*pi/2); % Transformed second initial guess
end

% Initial shape is ellipse due to initial guesses being negative
initial_guesses = [iv1,iv2];
initial_a = a_min./(1-initial_guesses.^2);
alpha = 2*acos( max(-1, min(1, initial_guesses)));
beta = long_way*2*asin(sqrt((s-c)/2./initial_a));

% Evaluate TOF with Lagrange
y12 = initial_a.*sqrt(initial_a).*((alpha - sin(alpha)) - (beta - sin(beta)) + 2*pi*nOrbits);

% Initial estimates for y
if nOrbits == 0
    y1 = log(y12(1)) - logTOF;
    y2 = log(y12(2)) - logTOF;
else
    y1 = y12(1) - TOF;
    y2 = y12(2) - TOF;
end

%% Solving for x and y
% This is the findxy function mentioned by Izzo in his paper.

% Newton's method
err = inf;
it = 0;
x_temp = 0;
while err > tolerance
    it = it + 1; % iterate
    x_temp = (x1*y2 - y1*x2)/(y2 - y1);

    if nOrbits == 0
        x = exp(x_temp) - 1;
    else
        x = atan(x_temp)*2/pi;
    end

    a = a_min/(1-x^2);

    % Handle ellipses and hyperbolas
    if x < 1 % Ellipse
        alpha = 2*acos( max(-1, min(1, x)));
        beta = long_way*2*asin(sqrt((s-c)/2/a));
    else % Hyperbola
        alpha = 2*acosh(x);
        beta = long_way*2*asinh(sqrt((s-c)/(-2*a)));
    end

    % Evaluate TOF with Lagrange
    if a > 0
        tof = a*sqrt(a)*((alpha-sin(alpha)) - (beta - sin(beta)) + 2*pi*nOrbits);
    else
        tof = -a*sqrt(-a)*((sinh(alpha) - alpha) - (sinh(beta) - beta));
    end

    % Update y
    if nOrbits == 0
        y_temp = log(tof) - logTOF;
    else
        y_temp = tof - TOF;
    end

    x1 = x2;
    x2 = x_temp;
    y1 = y2;
    y2 = y_temp;

    % Update error
    err = abs(x1 - x_temp);

    if it > 15
        check1 = true;
        break;
    end
end

% If Newton method fails, use Lancaster & Blanchard instead
if check1
    % ADD HERE
end

% Convert x
if nOrbits == 0
    x = exp(x_temp) - 1;
else
    x = atan(x_temp)*2*pi;
end

% Solve semimajor axis
a = a_min/(1-x^2);

% Solving psi and eta (see Izzo, p. 6)
if x < 1 % Ellipse
    alpha = 2*acos(max(-1, min(1, x2)));
    beta = long_way*2*asin(sqrt((s-c)/2/a));
    psi = (alpha-beta)/2;
    eta = sqrt(2*a*sin(psi)^2/s);
else % Hyperbola
    alpha = 2*acosh(x);
    beta = long_way*2*asinh(sqrt((c-s)/2/a));
        % Note: Earlier, this was (s-c) since there was an extra negative.
    psi = (alpha-beta)/2;
    eta = sqrt(-2*a*sinh(psi)^2/s);
end

% Compute v1_vec and v2_vec following remaining algorithm on Izzo p.15
ih = cr1r2_uv * long_way;
it1 = cross(ih,r1_vecND);
it2 = cross(ih,r2_vecUnit);
Vr1 = 1/(eta*sqrt(a_min))*(2*lambda*a_min - lambda - x*eta);
Vt1 = sqrt(r2_ND/(a_min*eta^2)*sin(dtheta/2)^2);
Vt2 = Vt1/r2_ND;
Vr2 = (Vt1 - Vt2)/tan(dtheta/2) - Vr1;
v1_vec = (Vr1*r1_vecND + Vt1*it1)*V;
v2_vec = (Vr2*r2_vecUnit + Vt2*it2)*V;


end


