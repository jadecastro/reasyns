function [Xdot] = CarKinematicsNL(t0,X)

global Xk Uk K t l Xss Uss Kss
persistent toff

n = 3;  % order of the Taylor expansion

x = X(1); y = X(2); th = X(3); phi = X(4);

% indx = find(min(abs(t - t0)) == abs(t - t0),1,'first');
% % weightedError = sum(repmat([1 1 0.1 0.01],size(Xk,1),1).*abs(Xk - repmat(X',size(Xk,1),1)),2);
% % indx = find(min(weightedError) == weightedError);
% K1 = K(:,:,indx);

% find the initial time offset
% indx = find(min(abs(t - t0)) == abs(t - t0));
if t0 == 0
    tmp = sum(repmat([1 1 0.001 0.001],size(Xk,1),1).*abs(Xk - repmat(X',size(Xk,1),1)),2);
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

% uV = -K1(1,:)*([x;y;th;phi] - Xk(indx,:)');
% uW = -K1(2,:)*([x;y;th;phi] - Xk(indx,:)');
uV = -4*K1(1,:)*([x;y;th;phi] - X1);
uW = -4*K1(2,:)*([x;y;th;phi] - X1);

% uV = Uk(indx,1) + uV;
% uW = Uk(indx,2) + uW;
% uV = max(min(uV,3),0);  % Limits angular velocity on (-1,1)
%uW = max(min(uW,10),-10);  % Limits angular velocity on (-1,1)

% Governing equations
% xdot1 = Uk(indx,1)*cos(Xk(indx,3));
% ydot1 = Uk(indx,1)*sin(Xk(indx,3));
% xdot1 = Uk(indx,1)*cos(th);
% ydot1 = Uk(indx,1)*sin(th);
% thdot1 = Uk(indx,1)*tan(phi)/l;
% phidot1 = Uk(indx,2);
xdot1 = U1(1)*cos(th);
ydot1 = U1(1)*sin(th);
thdot1 = U1(1)*tan(phi)/l;
phidot1 = U1(2);

xdot2 = uV*cosapprox(th,n);
ydot2 = uV*sinapprox(th,n);
% xdot2 = uV*cos(th);
% ydot2 = uV*sin(th);
thdot2 = uV*tan(phi)/l;
phidot2 = uW;

xdot = xdot1 + xdot2;
ydot = ydot1 + ydot2;
thdot = thdot1 + thdot2;
phidot = phidot1 + phidot2;

if t1 >= t(end)
    Xdot = [0; 0; 0; 0];
else
    Xdot = [xdot; ydot; thdot; phidot];
end
