function [Xdot] = CreateKinematicsCircleTverseCtrl(t,X)

global m J B r rw rho xd yd thd ud wd T1 T2 taun Bsa Pia Ka1 Ka2 cycle

global thd_tau

x = X(1); y = X(2); th = X(3);

thd_tau1 = thd_tau;

% define phase variable
tauna = atan2(y,x)/2/pi;
Bs_tau = interp1(taun,Bsa,tauna);
Pi_tau = interp1(taun,Pia,tauna);
K1_tau = interp1(taun,Ka1,tauna); K2_tau = interp1(taun,Ka2,tauna); 
xd_tau = interp1(taun,xd',tauna); yd_tau = interp1(taun,yd',tauna); 
thd_tau = interp1(taun,thd'-2*pi*cycle,tauna);
if abs(thd_tau - thd_tau1) > 2*pi*0.9
    cycle = cycle + 1;
    thd_tau = interp1(taun,thd'-2*pi*cycle,tauna);
end
    
Xd_tau = [xd_tau; yd_tau; thd_tau];

K_tau = [K1_tau; K2_tau];

% control law
ubar = -K_tau*[Pi_tau*([x;y] - Xd_tau(1:2)); (th - Xd_tau(3))];

ud_tau = interp1(taun,ud',tauna);
u_tau = [ud_tau; 0];
u_x = u_tau + ubar;


% Governing equations

xdot = u_x(1)*cos(th);
ydot = u_x(1)*sin(th);
thdot = u_x(2);

Xdot = [xdot; ydot; thdot];