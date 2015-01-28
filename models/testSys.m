function [Xdot] = testSys(t, X)

x1 = X(1);  x2 = X(2);

f = [(-x1^3 + x2);
    (-x2 - x1^2*x2)];

Xdot = f;
