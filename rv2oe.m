function coe = rv2oe(r_vec,v_vec,mu)

% Unit vectors
%I = [1;0;0]; 
%J = [0;1;0]; 
K = [0;0;1];

% Ensure that r_vec and v_vec are column vectors
[r_rows,r_cols] = size(r_vec);
[v_rows,v_cols] = size(v_vec);
if r_rows < r_cols
    r_vec = transpose(r_vec);
end
if v_rows < v_cols
    v_vec = transpose(v_vec);
end

r = norm(r_vec); %distance
v = norm(v_vec); %speed
v_r = dot(r_vec,v_vec)/r; %radial velocity
h_vec = cross(r_vec,v_vec); %specific angular momentum
h = norm(h_vec); %specific angular momentum magnitude
i = acos(h_vec(3)/h); %inclination in radians
N_vec = cross(K,h_vec); %node line
N = norm(N_vec);
if N_vec(2) >= 0
    Omega = acos(N_vec(1)/N); %RAAN
elseif N_vec(2) < 0
    Omega = 2*pi - acos(N_vec(1)/N); %RAAN
end
%RAAN = rad2deg(Omega); %RAAN in deg
e_vec = (1/mu)*((v^2-mu/r)*r_vec-r*v_r*v_vec); %eccentricity vector
e = norm(e_vec);
%e_vecT = e_vec';
if e_vec(3) >= 0
    omega = acos(dot(N_vec,e_vec)/(N*e)); %argument of periapse
elseif e_vec(3) < 0
    omega = 2*pi - acos(dot(N_vec,e_vec)/(N*e)); %argument of periapse
end
%AOP = rad2deg(omega); %AOP in deg

if e > 1.e-10
    TA = acos(dot(e_vec,r_vec)/e/r); %true anomaly
    if v_r < 0
        TA = 2*pi - TA;
    end
else
    cp = cross(N_vec,r_vec);
    if cp(3) >= 0
        TA = acos(dot(N_vec,r_vec)/N/r);
    else
        TA = 2*pi - acos(dot(N_vec,r_vec)/N/r);
    end
end
%TA = rad2deg(TA);
a = h^2/mu/(1-e^2);
%h_vecT = h_vec';
coe = [h e Omega i omega TA a];
end