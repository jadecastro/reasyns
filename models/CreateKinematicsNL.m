function [Xdot] = CreateKinematicsNL(t0,X)
% t0 is the sim time
% t is the vector of time indices for the base trajectory

global Xk Uk K t Xss Uss Kss
persistent toff

n = 3;  % order of the Taylor expansion

couplingCoeff = 0e-1;  % small coupling coefficient to make more amenable to sos computations

x = X(1); y = X(2); th = X(3);

% find the initial time offset
% indx = find(min(abs(t - t0)) == abs(t - t0));
if t0 == 0
    tmp = sum(repmat([1 1 0.00001],size(Xk,1),1).*abs(Xk - repmat(X',size(Xk,1),1)),2);
    indx = find(min(tmp) == tmp);
    toff = t(indx);
end
t1 = t0 + toff;

% K1 = K(:,:,indx);
K1 = Kss(t1);
X1 = Xss(t1);
U1 = Uss(t1);

% "unwrap" theta coordinate
% tmp = [abs(Xk(indx,3) - th); abs(Xk(indx,3) - (th + 2*pi)); abs(Xk(indx,3) - (th - 2*pi))];
tmp = [abs(X1(3) - th); abs(X1(3) - (th + 2*pi)); abs(X1(3) - (th - 2*pi))];
thetaCase = find(tmp == min(tmp),1,'first');
switch thetaCase
    case 2
        th = th + 2*pi;
    case 3
        th = th - 2*pi;
end

% uV = -1*K(1,:,indx)*([x;y;th] - Xk(indx,:)');
% uW = -1*K(2,:,indx)*([x;y;th] - Xk(indx,:)');
% uV = -100*K1(1,:)*([x;y;th] - X1);
% uW = -0.5*K1(2,:)*([x;y;th] - X1);
uV = [1 1 1].*K1(1,:)*([x;y;th] - X1);
uW = [1 1 1].*K1(2,:)*([x;y;th] - X1);

% uV = Uk(indx,1) + uV;
% uW = Uk(indx,2) + uW;
% uV = max(min(uV,3),0);  % Limits angular velocity on (-1,1)
%uW = max(min(uW,10),-10);  % Limits angular velocity on (-1,1)

% Governing equations
% xdot1 = Uk(indx,1)*cos(Xk(indx,3));
% ydot1 = Uk(indx,1)*sin(Xk(indx,3));
% xdot1 = Uk(indx,1)*cos(th);
% ydot1 = Uk(indx,1)*sin(th);
% thdot1 = Uk(indx,2);
xdot1 = U1(1)*cos(th);
ydot1 = U1(1)*sin(th);
thdot1 = U1(2);

% For debugging, use polynomial-approximated system - reflective of CreateKinematicsPoly
% xdot2 = uV*cosapprox(th,n);
% ydot2 = uV*sinapprox(th,n);
xdot2 = uV*cos(th);
ydot2 = uV*sin(th);
thdot2 = uW;

% Introduce a tiny amount of coupling - reflective of CreateKinematicsPoly
% xdot2 = xdot2 + couplingCoeff*x - couplingCoeff*y;
% ydot2 = ydot2 + couplingCoeff*y - couplingCoeff*x;
% thdot2 = thdot2 - couplingCoeff*th - couplingCoeff*x - couplingCoeff*y;

xdot = xdot1 + xdot2;
ydot = ydot1 + ydot2;
thdot = thdot1 + thdot2;

if t1 >= t(end)
    Xdot = [0; 0; 0];
else
    Xdot = [xdot; ydot; thdot];
end
