function chi = UniversalKepler(mu, DeltaT, r0, vr0, alpha)
% Solving the Kepler Equation using Newton's Method with Universal
% Variable, Chi
% Following Curtis algorithm

% chi, universal anomaly (e.g. km^0.5)
% mu, gravitational parameter (e.g. km^3/s^2)
% DeltaT, time since chi = 0 (e.g. s)
% r0, radial position (e.g. km) at chi = 0;
% vr0, radial velocity (e.g. km/s) when chi = 0;
% alpha, reciprocal of semimajor axis (e.g. 1/km)
% z, auxiliary variable (z = alpha*chi^2)
% C, value from Stumpff function
% S, value form Stumpff function
% n, number of iterations for convergence
% nMax, max number of iterations

% Requires stumpff.m

err = 1.e-8;
nMax = 1e5;

% Start value for chi
chi = sqrt(mu)*abs(alpha)*DeltaT;
%fprintf('Start value for chi = %f\n',chi);

% Iterate
n = 0;
ratio = 1;
while abs(ratio) > err && n <= nMax
    n = n + 1;
    z = alpha*chi^2;
    %fprintf('In Kepler: z = %f\n',z);
    [S,C] = stumpff(z);
    %fprintf('In Kepler: S = %f, C = %f\n',S,C);
    F = (r0*vr0)/sqrt(mu)*chi^2*C + (1-alpha*r0)*chi^3*S + r0*chi - sqrt(mu)*DeltaT;
    dFdchi = r0*vr0/sqrt(mu)*chi*(1-alpha*chi^2*S) + (1-alpha*r0)*chi^2*C + r0;
    %fprintf('F = %f, dFdchi = %f\n',F,dFdchi);

    ratio = F/dFdchi;
    chi = chi-ratio;
end

if n > nMax
    fprintf('Iterations = %i | F/dFdchi = %g\n',n, F/dFdchi);
end

%fprintf('End value for chi = %f\n',chi);
end
