function [S,C] = stumpff(z)
if z > 0
    S = (sqrt(z) - sin(sqrt(z)))/(sqrt(z^3));
    C = (1 - cos(sqrt(z)))/z;
elseif z < 0
    S = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt((-z)^3));
    C = (cosh(sqrt(-z)) - 1)/(-z);
else
    S = 1/6;
    C = 1/2;
end
end

%{
if z > 0
    S = (sqrt(z) - sin(sqrt(z)))/(sqrt(z))^3;
    C = (1 - cos(sqrt(z)))/z;
elseif z < 0
    S = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(-z))^3;
    C = (cosh(sqrt(-z)) - 1)/(-z);
else
    S = 1/6;
    C = 1/2;
end
%}




