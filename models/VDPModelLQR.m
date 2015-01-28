function [Xdot] = VDPModelLQR(t0,X)

global mu
global Xk Uk K t Xss Uss Kss
persistent toff

% chain of integrators

x = X(1); y = X(2);

if t0 == 0  % what point am I closest to at t=0?
    tmp = sum(repmat([1 1],size(Xk,1),1).*abs(Xk - repmat(X',size(Xk,1),1)),2);
    indx = find(min(tmp) == tmp);
    toff = t(indx);
end
t1 = t0 + toff;

% K1 = K(:,:,indx);
K1 = Kss(t1);
X1 = Xss(t1);
U1 = Uss(t1);

ux = -1*K1(1,:)*([x;y] - X1);
uy = -1*K1(2,:)*([x;y] - X1);

vx = U1(1);
vy = U1(2);

% Governing equations
xdot = y + vx + ux;
ydot = mu*(1 - x^2)*y - x + vy + uy;

if t1 >= t(end)
    Xdot = [0; 0];
else
    Xdot = [xdot; ydot];
end