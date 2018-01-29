function [t, p, pdot] = quinticSplines(q0,q1,q2,q3,v0,v3)


A = [ 1 0 0 0 0 0;
    1 1 1 1 1 1;
    1 2 4 8 16 32;
    1 3 9 27 81 243;
    0 1 0 0 0 0
    0 1 6 27 108 405];

B = [ q0; q1; q2; q3; v0; v3];

%solve for coeffs of polynomial
y = A\B

%reverse y to put in form for polyval
r=wrev(y)
t = 0:.1:3;

p = polyval(r,t);

k = polyder(r);
pdot = polyval(k,t);