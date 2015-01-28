function [Xdot] = CreateKinematicsNLAlt(t0,X)

global Xk Uk K t Xss Uss Kss
persistent toff

x = X(1); y = X(2); s = X(3); c = X(4);

% find the initial time offset
% indx = find(min(abs(t - t0)) == abs(t - t0));
if t0 == 0
    tmp = sum(repmat([1 1 0.001 0.001],size(Xk,1),1).*abs(Xk - repmat(X',size(Xk,1),1)),2);
    indx = find(min(tmp) == tmp,1,'first');
    toff = t(indx);
end
t1 = t0 + toff;

% K1 = K(:,:,indx);
K1 = Kss(t1);
X1 = Xss(t1);
U1 = Uss(t1);


% uV = -1*K(1,:,indx)*([x;y;th] - Xk(indx,:)');
% uW = -1*K(2,:,indx)*([x;y;th] - Xk(indx,:)');
% uV = -0.1*K1(1,:)*([x;y;th] - X1);
% uW = -0.1*K1(2,:)*([x;y;th] - X1);
uV = -1*K1(1,:)*([x;y;s;c] - X1);
uW = -1*K1(2,:)*([x;y;s;c] - X1);

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
xdot1 = U1(1)*c;
ydot1 = U1(1)*s;
sdot1 = U1(2)*c;
cdot1 = -U1(2)*s;

xdot2 = uV*c;
ydot2 = uV*s;
sdot2 = uW*c;
cdot2 = -uW*s;

xdot = xdot1 + xdot2;
ydot = ydot1 + ydot2;
sdot = sdot1 + sdot2;
cdot = sdot2 + cdot2;

if t1 >= t(end)
    Xdot = [0;0;0;0];
else
    Xdot = [xdot; ydot; sdot; cdot];
end
