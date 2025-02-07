% -----------------------
% Useage:
% computeCubicParameters( 0, 1, 0, 0, pi, 0)
% -----------------------
function x = computeCubicParameters( t0, tf, q0, v0, qf, vf)
y = [q0, v0, qf, vf]';
A = [1 t0 t0^2 t0^3; ...
0 1 2*t0 3*t0^2; ...
1 tf tf^2 tf^3;
0 1 2*tf 3*tf^2];
x = A\y;
end