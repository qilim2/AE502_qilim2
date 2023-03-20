function E = Laguerre(M,e)
% This function solves the Eccentric Anomaly given the Mean Anomaly and the
% eccentricity of the orbit using Laguerre's method.

% See AIAA-86-0084, An Improved Algorithm Due to Laguerre for the Solution
% of Kepler's Equation, Conway, 1986
% and also, Prussing and Conway Orbital Mechanics, p. 38

% Starting value
u = M + e;
n  = 5;
E  = (M*(1-sin(u)) + u*sin(M)) / (1 + sin(M) - sin(u));
for it1 = 1:5
    num = n*(E-e*sin(E));
    den = (1-e*cos(E)) + sign(1-e*cos(E))*(abs((n-1)^2*(1-e*cos(E))^2 ...
        - n*(n-1)*(E-e*sin(E))*(e*sin(E))))^0.5;
    E = num/den;
end
end